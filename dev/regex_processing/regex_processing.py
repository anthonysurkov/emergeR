import pandas as pd
from Bio import SeqIO
from itertools import islice
import gzip
import regex
import os
from tqdm import tqdm

infile = "r270x_preprocessed.fastq"
max_dist = 5

# template
tplt = (
    r"(CCTCTCCCAAGTCCACACAGAACGGGGCTGAAAGCCGGGTTTCTTTCTTTC"
    r"CCGNNNNNNNNNCCCGTTTGCCCGTAGAGTCGCTGTTC)"
)
tplt = tplt.replace("U", "T").replace("A", "[AG]").replace("N", "[ACGT]")
tplt = tplt + "{e<=" + str(max_dist) + "}"
tplt = regex.compile(tplt, regex.BESTMATCH)

# target and variable capture
trgt = regex.compile(r"(?:GCTG){e<=1}([AG]{3})(?:GCCG){e<=1}")
varb = regex.compile(r"(?:TTTCTTTCTTTC){e<=1}.*?([ACGT]{9})(?:CCCGTTTGCCCG){e<=2}")

def smart_open(path):
    with open(path, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":  # gzip magic bytes
        return gzip.open(path, "rt", encoding="latin-1")
    return open(path, "r", encoding="latin-1")


def chunk_fastq(path, chunk_size=5000):
    total_size = os.path.getsize(path)
    with smart_open(path) as f, tqdm(
        total=total_size, unit="B", unit_scale=True, desc="Reading FASTQ"
    ) as pbar:
        records = []
        while True:
            start_pos = f.tell()
            # Each read = 4 lines
            for _ in range(chunk_size):
                id_line = f.readline()
                if not id_line:
                    break
                seq = f.readline().strip()
                _ = f.readline()
                qual = f.readline().strip()
                records.append((id_line[1:].strip(), seq, qual))
            if not records:
                break
            yield pd.DataFrame(records, columns=["id", "seq", "qual"])
            records = []
            pbar.update(f.tell() - start_pos)

def regex_match(seq):
    m = tplt.search(seq)
    if not m:
        return None, None
    read = m.group()
    target_m = trgt.search(read)
    if not target_m:
        return None, None
    n10_m = varb.search(read)
    if not n10_m:
        return None, None
    return n10_m.group(1), target_m.group(1)

def casey_match(seq, template, expected, flanks, n10_region, edit_region):
    def subset(seq, start, end):
        return str(seq[start-1:end])

    concat = (
        subset(seq, flanks[0][0], flanks[0][1]) +
        subset(seq, flanks[1][0], flanks[1][1]) +
        subset(seq, flanks[2][0], flanks[2][1]) +
        subset(seq, flanks[3][0], flanks[3][1])
    )
    if concat != expected:
        return None, None

    n10 = subset(seq, n10_region[0], n10_region[1])
    target = subset(seq, edit_region[0], edit_region[1])

    return n10, target

records = []
for df in chunk_fastq(infile, chunk_size=10000):
    for s in df["seq"]:
        n10, target = regex_match(s)
        if n10 and target:
            records.append((n10, target))

if not records:
    print("No matches found")
else:
    table = pd.DataFrame(records, columns=["n10", "target"])
    ctab = pd.crosstab(table["n10"], table["target"])
    tot = int(ctab.values.sum())

    #print(ctab.to_string())
    print(f"\nNumber of unique n10s: {ctab.shape[0]}")
    print(f"Total reads matched: {tot}")

    ctab.to_csv("regex_out.csv")
    print("Saved to regex_out.csv")

