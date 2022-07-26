import pandas as pd
import glob
import os
def assemble_raw_counts(path="."):
    file = "gene_count.txt"
    count_files = []

    for root,dirs,files in os.walk(path):
        if file in files:
            count_files.append(root + "/" + file)
    count_files = sorted(count_files)
    for f, file in enumerate(count_files) :
        df = pd.read_table(file, comment = "#")
        name = [df.columns[len(df.columns)-1].split("/")][0][1]
        if f == 0:
            final_df = df.iloc[:, [0,5,7,8,9]]
            final_df.rename(columns={df.columns[len(df.columns)-1]: name}, inplace=True)
        else:
            final_df[name] = df.iloc[:, 9]
    return final_df

df = assemble_raw_counts(snakemake.params.countDir)
df.to_csv(snakemake.params.countDir + "/count-tables/raw_counts.tab",sep="\t", index=False)

