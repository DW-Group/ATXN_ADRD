import pandas as pd
import gwaslab as gl
import os
import matplotlib as mpl

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

# -----------------------------
# input / output
# -----------------------------

infile = "results/combined_mht_sumstats_top1_per_snp.txt"
outfile = "figures/fig5a_manhattan.pdf"

# -----------------------------
# load
# -----------------------------

df = pd.read_csv(infile, sep="\t")

df["Chr"] = pd.to_numeric(df["Chr"], errors="coerce")
df["Pos"] = pd.to_numeric(df["Pos"], errors="coerce")
df["Pval"] = pd.to_numeric(df["Pval"], errors="coerce")

df = df.dropna(subset=["Chr", "Pos", "Pval"]).copy()
df = df[(df["Pval"] > 0) & (df["Pval"] <= 1)].copy()

df["Chr"] = df["Chr"].astype(int)
df["Pos"] = df["Pos"].astype(int)

df["SNPID"] = df["Chr"].astype(str) + ":" + df["Pos"].astype(str)

df["rsid"] = df["rsid"].fillna(".").astype(str)
df["Trait"] = df["Trait"].fillna("NA").astype(str)

df["LABEL"] = df["Trait"] + " | " + df["rsid"]

df = df.sort_values("Pval", ascending=True).copy()

# -----------------------------
# annotate genome-wide significant hits
# -----------------------------

sig_df = (
    df.loc[df["Pval"] < 5e-8]
      .copy()
)

sig_df = sig_df.sort_values(["Trait", "Chr", "Pos", "Pval"]).copy()

# define loci within each protein
locus_window = 500_000

def assign_locus(group):
    group = group.sort_values(["Chr", "Pos", "Pval"]).copy()
    locus_ids = []
    current_locus = 0
    prev_chr = None
    prev_pos = None
    
    for _, row in group.iterrows():
        if prev_chr is None:
            current_locus += 1
        elif row["Chr"] != prev_chr:
            current_locus += 1
        elif row["Pos"] - prev_pos > locus_window:
            current_locus += 1
        
        locus_ids.append(current_locus)
        prev_chr = row["Chr"]
        prev_pos = row["Pos"]
    
    group["locus_id"] = locus_ids
    return group

sig_df = (
    sig_df.groupby("Trait", group_keys=False)
          .apply(assign_locus)
)

# keep only top SNP per protein per locus
anno_df = (
    sig_df.sort_values(["Trait", "locus_id", "Pval"], ascending=[True, True, True])
          .drop_duplicates(subset=["Trait", "locus_id"], keep="first")
          .copy()
)

anno_set = anno_df["SNPID"].tolist()
anno_alias = dict(zip(anno_df["SNPID"], anno_df["LABEL"]))

# -----------------------------
# load into gwaslab
# -----------------------------

mysumstats = gl.Sumstats(
    df,
    snpid="SNPID",
    chrom="Chr",
    pos="Pos",
    p="Pval",
    build="19",
    verbose=False
)

mysumstats.fix_chr(verbose=False)

# -----------------------------
# plot
# -----------------------------

mysumstats.plot_mqq(
    mode="m",
    skip=0,
    cut=None,
    anno=True,
    anno_set=anno_set,
    anno_alias=anno_alias,
    anno_line_color="grey",
    anno_line_style="--",
    anno_sig_level=5e-8,
    colors = ["#9FB2C6", "#5C779C"],
    scatter_kwargs={"s": 4, "alpha": 0.9, "edgecolor": "none"},
    sig_line_color="black",
    sig_line_width=0.2,
    check=False,
    verbose=False,
    figsize=(14, 4),
    dpi=1000,
    save=outfile
)

print(f"Saved plot to: {outfile}")