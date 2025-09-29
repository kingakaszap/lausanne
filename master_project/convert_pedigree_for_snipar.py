import pandas as pd

infile = "data/coreccted_consensus_pedigree.txt"
outfile = "data/snipar_pedigree.txt"

df = pd.read_csv(infile, sep=",", dtype=str)

df =df.rename(columns={
    "RingId": "IID",
    "MumId": "MOTHER_ID",
    "DadId": "FATHER_ID"
})

for col in ["MOTHER_ID", "FATHER_ID"]:
    df[col] = df[col].fillna("0").replace({"NA":"0", "Na":"0", "na":"0", "":"0"})

def make_fid(row):
    if row["MOTHER_ID"] != "0":
        return row["MOTHER_ID"]
    if row["FATHER_ID"] !="0":
        return row["FATHER_ID"]
    return row["IID"]

df["FID"] = df.apply(make_fid, axis=1)

out = df[["FID", "IID", "FATHER_ID", "MOTHER_ID"]]

out.to_csv(outfile, sep=" ", index=False)