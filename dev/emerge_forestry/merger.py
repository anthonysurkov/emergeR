import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import re
from networkx.drawing.nx_agraph import graphviz_layout
from collections import Counter, deque, defaultdict

# --- Helpers --- #
def parse_token(seq):
    matches = re.findall(r'([ACGU])(\d+)', seq)
    return {int(pos): base for base, pos in matches}

def build_kmer_index(seq_pop, seq_col="n10", edit_col="map", kmax=6):
    index = {}
    for _, row in seq_pop.iterrows():
        seq = row[seq_col].replace("T", "U")
        edit = row[edit_col]

        # tokenize with position info
        tokens = [f"{c}{i}" for i, c in enumerate(seq)]
        for k in range(2, kmax+1):
            for i in range(len(tokens)-k+1):
                kmer = "".join(tokens[i:i+k])
                positions = tuple(int(p) for _, p in re.findall(r'([ACGU])(\d+)', kmer))
                motif = "".join(b for b, _ in re.findall(r'([ACGU])(\d+)', kmer))
                key = (positions, motif)
                index.setdefault(key, []).append(edit)
    # average
    return {k: sum(v)/len(v) for k, v in index.items()}

def find_avg_edit(node, kmer_index):
    parsed = parse_token(node.get_seq())   # {6:"G", 7:"A", 8:"U"}
    if not parsed:
        return float("nan")
    positions = tuple(sorted(parsed.keys()))
    motif = "".join(parsed[p] for p in positions)
    return kmer_index.get((positions, motif), float("nan"))

def append_avg_edits_to_tree(node, seq_index):
    node.set_editing(find_avg_edit(node, seq_index))
    for child in node.get_children():
        append_avg_edits_to_tree(child, seq_index)

def append_avg_edits_to_forest(forest, seq_index):
    for root in forest:
        append_avg_edits_to_tree(root, seq_index)

def convert_to_rna(dna):
    rna = []
    for char in dna:
        if char == "T":
            rna.append("U")
        else:
            rna.append(char)
    rna = "".join(rna)
    return rna

# --- Visualizers --- #
def print_tree(node, depth=0):
    id = node.get_id()
    seq = node.get_seq()
    length = node.get_len()
    rank = node.get_rank()
    prevalence = node.get_prevalence()
    editing = node.get_editing()
    if rank:
        print("  " * depth +
             f"Kmer {id}, len {length}, {seq}, rank #{rank} "
             f"(per-{length}-mer), count {prevalence}, "
             f"avg editing {editing}"
        )
    else:
        print("  " * depth +
             f"Kmer {id}, len {length}, {seq}, "
             f"(per-{length}-mer), count {prevalence}, "
             f"avg editing {editing}"
        )
    for child in node.get_children():
        print_tree(child, depth+1)

def forest_to_nx(forest):
    G = nx.DiGraph()
    for root in forest:
        stack = [root]
        while stack:
            node = stack.pop()
            for child in node.get_children():
                parent_label = f"{node.get_seq()}"
                child_label = f"{child.get_seq()}"
                G.add_edge(child_label, parent_label)  # small -> big

                # add attributes for coloring/size
                G.nodes[parent_label]["rank"] = node.get_rank()
                G.nodes[parent_label]["count"] = node.get_prevalence()
                G.nodes[parent_label]["editing"] = node.get_editing()

                G.nodes[child_label]["rank"] = child.get_rank()
                G.nodes[child_label]["count"] = child.get_prevalence()
                G.nodes[child_label]["editing"] = child.get_editing()

                stack.append(child)
    return G

# --- Node definition --- #
class Node:
    def __init__(self, id, seq=None, rank=None, parent=None, prevalence=0):
        self.id = id
        self.seq = seq
        self.rank = rank
        self.parents = []
        self.children = []
        self.prevalence = prevalence
        self.editing = 0
        if parent:
            self.parents.append(parent)

    def set_child(self, child):
        self.children.append(child)
    def set_parent(self, parent):
        self.parents.append(parent)
    def set_seq(self, seq):
        self.seq = seq
    def set_id(self, id):
        self.id = id
    def set_rank(self, rank):
        self.rank = rank
    def set_prevalence(self, prevalence):
        self.prevalence = prevalence
    def set_editing(self, editing):
        self.editing = editing

    def get_children(self):
        return self.children
    def get_parents(self):
        return self.parents
    def get_seq(self):
        return self.seq
    def get_len(self):
        if not self.seq:
            return 0
        return int(len(self.seq) / 2)
    def get_id(self):
        return self.id
    def get_rank(self):
        return self.rank
    def get_prevalence(self):
        return self.prevalence
    def get_editing(self):
        return self.editing

# --- BPE implementation --- #
class EmergeTokenizer:
    def __init__(self, bases, seqs, kmax):
        if not bases:
            raise ValueError(
                "Please provide valid RNA bases (non-canonical OK)"
            )
        if not seqs.any():
            raise ValueError(
                "Please provide a population of sequences as a list"
                " (alignment not necessary, must be the same length)"
            )
        if not kmax or kmax < 2:
            raise ValueError(
                "Kmax should be greater than or equal to 2"
            )

        #same_length = all(len(s) == len(seqs[0]) for s in seqs)
        #if not same_length:
        #    raise ValueError("Not all sequences same length")


        self.rna_chars = self.define_rna_chars(bases)
        self.vocab = {i: char for i, char in enumerate(self.rna_chars)}
        self.inverse_vocab = {char: i for i, char in self.vocab.items()}
        self.kmax = kmax
        self.seqs = seqs

        self.merges = {}
        self.global_ranks = {}
        self.kmer_ranks = {}
        self.kmer_bins = defaultdict(list)
        self.counts = Counter()

    def define_rna_chars(self, bases):
        rna_chars = [f"{base}{i}" for base in bases for i in range(0, 10)]
        return rna_chars

    def encode(self, seqs=None, vocab_size=1000, special="E", to_rna=True):
        if not seqs:
            if not self.seqs.any():
                raise ValueError("Please provide sequences")
            seqs = self.seqs

        processed_rna = []
        for seq in seqs:
            if to_rna:
                seq = convert_to_rna(seq)
            for i, char in enumerate(seq):
                processed_rna.append(f"{char}{i}")
            processed_rna.append("E")

        if special:
            for token in special:
                if token not in self.inverse_vocab:
                    new_id = len(self.vocab)
                    self.vocab[new_id] = token
                    self.inverse_vocab[token] = new_id

        token_ids = [self.inverse_vocab[token] for token in processed_rna]

        global_rank = 0
        for new_id in range(len(self.vocab), vocab_size):
            pair_id = self.find_freq_pair(token_ids, "most", self.vocab)
            if pair_id is None:
                break

            p0, p1 = pair_id
            merged_token = self.vocab[p0] + self.vocab[p1]
            if len(merged_token) > 2 * self.kmax:
                continue

            token_ids = self.replace_pair(token_ids, pair_id, new_id)
            global_rank += 1
            self.merges[pair_id] = new_id
            self.vocab[new_id] = merged_token
            self.inverse_vocab[merged_token] = new_id
            self.global_ranks[(p0, p1)] = global_rank

            k = len(merged_token)
            self.kmer_bins[k].append(new_id)

        for k, ids in self.kmer_bins.items():
            for i, tid in enumerate(ids, start=1):
                self.kmer_ranks[tid] = i

        self.counts = self.count_kmers(seqs, self.kmax)

    def build_forest(self):
        nodes = {}

        def get_node(id):
            if id not in nodes:
                nodes[id] = Node(id=id, seq=str(id))
            return nodes[id]

        for (left_id, right_id), parent_id in self.merges.items():
            parent = get_node(parent_id)
            left = get_node(left_id)
            right = get_node(right_id)

            parent.set_child(left)
            parent.set_child(right)
            left.set_parent(parent)
            right.set_parent(parent)

            left.set_seq(self.vocab[left_id])
            right.set_seq(self.vocab[right_id])
            parent.set_seq(self.vocab[parent_id])

            parent.set_rank(self.kmer_ranks.get(parent_id))

            parent.set_prevalence(self.counts.get(parent.seq, 0))
            left.set_prevalence(self.counts.get(left.seq, 0))
            right.set_prevalence(self.counts.get(right.seq, 0))

        all_children = {c for (l, r) in self.merges for c in (l, r)}
        roots = [get_node(pid) for (l, r), pid in self.merges.items()
                 if pid not in all_children]
        return roots

    def count_kmers(self, seqs, kmax):
        counts = Counter()
        for seq in seqs:
            seq = convert_to_rna(seq)
            tokens = [f"{char}{i}" for i, char in enumerate(seq)]
            for k in range(2, kmax + 1):
                for i in range(len(tokens) - k + 1):
                    kmer = "".join(tokens[i:i+k])
                    counts[kmer] += 1
        return counts

    def find_freq_pair(self, token_ids, mode="most", vocab=None):
        pairs = Counter(zip(token_ids, token_ids[1:]))

        # TODO: replace "E" hardcoding with all special characters
        if vocab is not None:
            pairs = {p: c for p, c in pairs.items()
                     if "E" not in vocab[p[0]] and "E" not in vocab[p[1]]}
        if mode == "most":
            return max(pairs.items(), key=lambda x: x[1])[0]
        elif mode == "least":
            return min(pairs.items(), key=lambda x: x[1])[0]
        else:
            raise ValueError("Invalid mode. Choose 'most' or 'least'.")

    def replace_pair(self, token_ids, pair_id, new_id):
        dq = deque(token_ids)
        replaced = []
        while dq:
            current = dq.popleft()
            if dq and (current, dq[0]) == pair_id:
                replaced.append(new_id)
                dq.popleft()
            else:
                replaced.append(current)
        return replaced

def temp(top_seqs, all_data, the_filename, title, edit_col):
    bases = ["A", "G", "C", "U"]
    Kmax = 6

    tokenizer1 = EmergeTokenizer(bases, top_seqs, Kmax)
    tokenizer1.encode()
    kmer_forest = tokenizer1.build_forest()

    seq_index = build_kmer_index(all_data, seq_col="n10", edit_col=edit_col)
    append_avg_edits_to_forest(kmer_forest, seq_index)

    for root in kmer_forest:
        print_tree(root)

    """
    # visualization example: color nodes by count, size by count
    G = forest_to_nx(kmer_forest)

    edits = nx.get_node_attributes(G, "editing")
    node_colors = [edits.get(n, 0) for n in G.nodes()]
    node_sizes = [100 + 5*nx.get_node_attributes(G, "count").get(n, 0) for n in G.nodes()]

    pos = graphviz_layout(G, prog="dot")

    plt.figure(figsize=(12, 8))
    ax = plt.gca()   # <-- explicitly grab current axes

    nodes = nx.draw(
        G, pos,
        with_labels=True,
        node_color=node_colors,
        node_size=node_sizes,
        cmap=plt.cm.Reds,
        vmin=min(node_colors),   # scale dynamically
        vmax=max(node_colors),
        font_size=8,
        arrows=True,
        ax=ax
    )

    sm = plt.cm.ScalarMappable(cmap=plt.cm.Reds,
                               norm=plt.Normalize(vmin=min(node_colors),
                                                  vmax=max(node_colors)))
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label="Avg editing level")

    plt.show()
    """

    from pyvis.network import Network
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors

    def forest_to_pyvis(forest, filename="forest.html"):
        # make interactive directed network
        net = Network(height="750px", width="100%", directed=True, notebook=False)
        net.toggle_physics(False)         # keep Graphviz layout
        net.set_edge_smooth('straight')   # no spaghetti

        # get layout from graphviz
        import networkx as nx
        from networkx.drawing.nx_agraph import graphviz_layout

        G = nx.DiGraph()
        for root in forest:
            stack = [root]
            while stack:
                node = stack.pop()
                for child in node.get_children():
                    G.add_edge(child.get_seq(), node.get_seq())
                    stack.append(child)

        pos = graphviz_layout(G, prog="dot")

        # color mapping: gray→red by editing
        edit_values = []
        for root in forest:
            stack = [root]
            while stack:
                node = stack.pop()
                if node.get_editing() is not None:
                    edit_values.append(node.get_editing())
                stack.extend(node.get_children())

        if edit_values:
            norm = mcolors.Normalize(vmin=min(edit_values), vmax=max(edit_values))
            cmap = cm.get_cmap("Reds")

        # add nodes + edges with attributes
        for n in G.nodes():
            # find node object back
            # (assuming seqs are unique)
            edit = None
            prevalence = None
            rank = None
            for root in forest:
                stack = [root]
                while stack:
                    node = stack.pop()
                    if node.get_seq() == n:
                        edit = node.get_editing()
                        prevalence = node.get_prevalence()
                        rank = node.get_rank()
                    stack.extend(node.get_children())

            color = "#cccccc"
            if edit is not None:
                rgba = cmap(norm(edit))
                color = mcolors.to_hex(rgba)

            title = f"""Seq: {n}
            Rank: {rank}
            Count: {prevalence}
            Avg editing: {edit:.3f}
            """

            net.add_node(
                n,
                label=n,
                title=title,
                color=color,
                size=10 + (prevalence or 1) * 0.1,
            )

        for u, v in G.edges():
            net.add_edge(u, v)

        # lock positions from graphviz
        for n, (x, y) in pos.items():
            net.nodes[list(G.nodes()).index(n)]['x'] = x
            net.nodes[list(G.nodes()).index(n)]['y'] = -y  # flip Y for pyvis

        net.write_html(filename, open_browser=True)

    def add_titles_to_html(html_file, title="Graph Title", subtitle="Graph Subtitle"):
        with open(html_file, "r") as f:
            html = f.read()

        insertion = f"""
        <div style="text-align:center; margin-bottom: 20px;">
            <h2 style="margin:0;">{title}</h2>
            <p style="margin:0; color:dark-gray; font-size:14px;">{subtitle}</p>
        </div>
        """

        html = html.replace("<body>", f"<body>{insertion}", 1)

        with open(html_file, "w") as f:
            f.write(html)

    forest_to_pyvis(kmer_forest, the_filename)
    subtitle = """
        0-indexed. NNNNNNNNNN = 0123456789<br>
        e.g. G6A7U8 is NNNNNNGAUN<br>
        e.g. U0U1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;is UUNNNNNNNN<br><br>
        Mouse over each node for more information.<br>
        Rank refers to the order in which byte-pair encoding (BPE) identified
         the motifs, adjusted for k-mer size. For example, rank 1 means that
         it was the first k-mer of its length to be identified, making it the
         most frequent motif of its length.<br>
        Count is the number of times the motif appears in the high-editing
         population. The size of a node also roughly corresponds to
         its prevalence.<br>
        Average editing is taken across the entire population of sequences,
         not just the high-editing population.
    """
    add_titles_to_html(the_filename, title, subtitle)

def main():
    edit_col = "mle"

    df = pd.read_csv("data/R270X_stats.csv")
    df = df.rename(columns={"N10": "n10"})
    df = df[["n10", edit_col]]

    df = df.sort_values(by=edit_col, ascending=False)
    df_top = df[df[edit_col] > 0.90]
    top_seqs = df_top["n10"]
    print(df_top.shape[0])

    outfile = "forest_r270x_map_over_060.html"
    title = "R270X, MAP > 0.60 (n = 1222)"
    #temp(top_seqs, df, outfile, title, edit_col)

if __name__ == "__main__":
    main()
