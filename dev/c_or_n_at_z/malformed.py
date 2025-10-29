import pandas as pd
import regex

df = pd.read_csv("r270x_cstatus.csv")

template = r"([AGCT]{4})(?:[AGCT])([AGCT]{4})"
template = regex.compile(template)

n8_dict = {}
for seq in df['n10']:
    m = template.match(seq)
    if not m:
        continue
    key = m.group(1) + m.group(2)
    if key not in n8_dict:
        n8_dict[key] = 1
    else:
        n8_dict[key] += 1

problems = {k: v for k, v in n8_dict.items() if v > 1}
problem_len = len(problems)

# problem = my affectionate name for "an N9 region which has the same N1-N4 and
#           N6-N9 regions. This would evidence the middle base being read
#           differently
#
# 0 problems detected
print(f"Total R270X length: {df.shape[0]}")
print(f"We have {problem_len} total problems in R270X")

