# input data (copy-paste your table here)
data = """
1	340.030995	1.000	1.000
2	171.147568	1.987	0.993
4	89.954681	3.780	0.945
8	44.429740	7.653	0.957
16	23.215947	14.646	0.915
32	15.129021	22.475	0.702
"""

rows = [line.split() for line in data.strip().splitlines()]

# print markdown table
print("| #processes |    Time     | Speedup | Efficiency |")
print("|------------|-------------|---------|------------|")

for p, t, s, e in rows:
    print(f"| {p:<10} | {float(t):<11.6f} | {float(s):<7.3f} | {float(e):<10.3f} |")
