# Define an isomorphism as a set of sequence identities such that each carries
# an intentionally-fixed base, misread in deamination post-processing.
#
# e.g. for NNNNZNNNN, we read AAAACTTTT and AAAAGTTTT. These are isomorphic, as
# C and G were Z during the reaction.
# We pool the associated reads of isomorphic identities to estimate true
# editing.

import pandas as pd
import regex
from tqdm import tqdm

# make this gooder later
df = pd.read_csv("../../data/r270x_regexed.csv")
df['k'] = df['GAA']
df['n'] = df['AAA'] + df['GAA']
df = df[['n10', 'n', 'k']]
df['mle'] = df['k'] / df['n']

template = regex.compile(r"([AGCT]{4})(?:[AGCT])([AGCT]{4})")

n8_dict = {}
for seq in tqdm(df['n10'], desc="Pooling isomorphs"):
    m = template.match(seq)
    if not m:
        continue
    key = m.group(1) + 'Z' + m.group(2)
    if key not in n8_dict:
        n8_dict[key] = {
            "n": df.loc[df['n10'] == seq, 'n'].iloc[0],
            "k": df.loc[df['n10'] == seq, 'k'].iloc[0]
        }
    else:
        n8_dict[key]['n'] += df.loc[df['n10'] == seq, 'n'].iloc[0]
        n8_dict[key]['k'] += df.loc[df['n10'] == seq, 'k'].iloc[0]

df = (
    pd.DataFrame.from_dict(n8_dict, orient='index', columns=['n', 'k'])
        .rename_axis('n10')
        .reset_index()
)

df.to_csv('r270x_no_isomorphs.csv')
