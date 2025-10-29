import pandas as pd

infile_regex = "benchmark_regex.csv"
infile_casey1 = "benchmark1_casey.csv"
infile_casey2 = "benchmark2_casey.csv"
infile_casey3 = "benchmark3_casey.csv"

df_rx = pd.read_csv(infile_regex)
df_c1 = pd.read_csv(infile_casey1)
df_c2 = pd.read_csv(infile_casey2)
df_c3 = pd.read_csv(infile_casey3)
df_c1 = df_c1.rename(columns={df_c1.columns[0]: 'n10'})
df_c2 = df_c2.rename(columns={df_c2.columns[0]: 'n10'})
df_c3 = df_c3.rename(columns={df_c3.columns[0]: 'n10'})

# n10 = number of unique n10s
rx_n10 = df_rx['n10'].shape[0]
c1_n10 = df_c1['n10'].shape[0]
c2_n10 = df_c2['n10'].shape[0]
c3_n10 = df_c3['n10'].shape[0]

# ont = on-target edit reads
rx_ont = df_rx['AAA'].sum() + df_rx['GAA'].sum()
c1_ont = df_c1['AAA'].sum() + df_c1['GAA'].sum()
c2_ont = df_c2['AAA'].sum() + df_c2['GAA'].sum()
c3_ont = df_c3['AAA'].sum() + df_c3['GAA'].sum()

# tot = total reads
rx_tot = df_rx.select_dtypes(include='number').sum().sum()
c1_tot = df_c1.select_dtypes(include='number').sum().sum()
c2_tot = df_c2.select_dtypes(include='number').sum().sum()
c3_tot = df_c3.select_dtypes(include='number').sum().sum()

# chg = change
c1_n10_chg = ((rx_n10 - c1_n10) / rx_n10) * 100
c2_n10_chg = ((rx_n10 - c2_n10) / rx_n10) * 100
c3_n10_chg = ((rx_n10 - c3_n10) / rx_n10) * 100
c1_ont_chg = ((rx_ont - c1_ont) / rx_ont) * 100
c2_ont_chg = ((rx_ont - c2_ont) / rx_ont) * 100
c3_ont_chg = ((rx_ont - c3_ont) / rx_ont) * 100
c1_tot_chg = ((rx_tot - c1_tot) / rx_tot) * 100
c2_tot_chg = ((rx_tot - c2_tot) / rx_tot) * 100
c3_tot_chg = ((rx_tot - c3_tot) / rx_tot) * 100

print("\nBenchmark on first 1000 R270X rawdata lines\n")
print("Regex:")
print(f"\tUnique N10s: {rx_n10}")
print(f"\tTotal reads: {rx_tot}")
print(f"\tAAA+GAA reads: {rx_ont}")
print("Flank Set 1:")
print(f"\tUnique N10s: {c1_n10} ({c1_n10_chg:.2f}% fewer than regex)")
print(f"\tTotal reads: {c1_tot} ({c1_tot_chg:.2f}% fewer)")
print(f"\tAAA+GAA reads: {c1_ont} ({c1_ont_chg:.2f}% fewer)")
print("Flank Set 2:")
print(f"\tUnique N10s: {c2_n10} ({c2_n10_chg:.2f}% fewer)")
print(f"\tTotal reads: {c2_tot} ({c2_tot_chg:.2f}% fewer)")
print(f"\tAAA+GAA reads: {c2_ont} ({c2_ont_chg:.2f}% fewer)")
print("Flank Set 3:")
print(f"\tUnique N10s: {c3_n10} ({c3_n10_chg:.2f}% fewer)")
print(f"\tTotal reads: {c3_tot} ({c3_tot_chg:.2f}% fewer)")
print(f"\tAAA+GAA reads: {c3_ont} ({c3_ont_chg:.2f}% fewer)")
print()

