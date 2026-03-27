import pandas as pd

df = pd.read_csv('matrix_DPP4.csv', sep='\t')

print(df.shape)
print(df.head())

# Each residue appears once per interaction type (8 interaction types total).
# With 7 residues, the columns for a residue are spaced by +7 each block.

df['Glu205'] = pd.DataFrame(df.iloc[:, [1, 8, 15, 22, 29, 36, 43, 50]]).apply(any, axis=1)
df['Glu206'] = pd.DataFrame(df.iloc[:, [2, 9, 16, 23, 30, 37, 44, 51]]).apply(any, axis=1)
df['Ser209'] = pd.DataFrame(df.iloc[:, [3, 10, 17, 24, 31, 38, 45, 52]]).apply(any, axis=1)
df['Phe357'] = pd.DataFrame(df.iloc[:, [4, 11, 18, 25, 32, 39, 46, 53]]).apply(any, axis=1)
df['Arg358'] = pd.DataFrame(df.iloc[:, [5, 12, 19, 26, 33, 40, 47, 54]]).apply(any, axis=1)
df['Tyr547'] = pd.DataFrame(df.iloc[:, [6, 13, 20, 27, 34, 41, 48, 55]]).apply(any, axis=1)
df['Tyr662'] = pd.DataFrame(df.iloc[:, [7, 14, 21, 28, 35, 42, 49, 56]]).apply(any, axis=1)

# Count how many of the 7 binding-site residues are hit (True) for each ligand
df['AllInteractions'] = (
    pd.Series(df[['Glu205','Glu206','Ser209','Phe357','Arg358','Tyr547','Tyr662']].values.tolist())
      .apply(lambda x: sum(x))
)

df.to_csv('matrix_DPP4_interactions_count.csv', sep='\t', index=False)

