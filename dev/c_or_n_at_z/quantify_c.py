import pandas as pd
from collections import Counter
import regex

infile = "../../data/r270x_regexed.csv"
df = pd.read_csv(infile)
#df_high = df[(df['mle'] > 0.30)]
c_regex = regex.compile(r"(?:[ACTG]){4}(C)(?:[ACTG]){4}")

def get_c_status(x):
    if bool(c_regex.search(str(x))):
        return True
    return False

df["c_status"] = df["n10"].apply(get_c_status)
#df_high["c_status"] = df["n10"].apply(get_c_status)

# for all
count = df["c_status"].sum()
total = df["c_status"].shape[0]
percent = (count / total) * 100
print(f"{percent}% of Z position is C (all N9s, total n = {total})")

# for high
#count_h = df_high["c_status"].sum()
#total_h = df_high["c_status"].shape[0]
#percent_h = (count_h / total_h) * 100
#print(f"{percent_h}% of Z position is C (mle > 0.30, total n = {total_h})")

#df = df[["n10", "mle", "map", "c_status"]]
#df.to_csv("r270x_cstatus.csv")
