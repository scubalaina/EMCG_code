import pandas as pd
import os, sys, re
from Bio import SeqIO

# run by command:
# python 3_parse_virdetect.py

sag_lens = pd.read_csv("test_untrimmed_contig_lens.tsv",sep="\t",header=0,index_col=False,low_memory=False)

sag_lens["genome"] = sag_lens["contig"].str.split("_").str[0:-1]
sag_lens["genome"] = sag_lens["genome"].str.join("_")
sag_lens = sag_lens[["genome","contig","length"]]

vs2_files = []
for r, d, f in os.walk("test_vs2_out/"):
    for filename in f:
        if filename.endswith("_vs2-final-viral-score.tsv"):
            fullfile = os.path.join(r,filename)
            vs2_files.append(fullfile)

vs2df = pd.concat([pd.read_csv(fp, sep="\t",low_memory=False) for fp in vs2_files])
vs2df["contig"] = vs2df["seqname"].str.split("\|\|").str[0]
vs2df.drop(["seqname"],inplace=True,axis=1)

newcol_vs2 = []
for x in vs2df.columns:
    if x == "contig":
        newcol_vs2.append(x)
    else:
        new = "vs2_" + x
        newcol_vs2.append(new)

vs2df.columns = newcol_vs2

dvf1_files = []
for r, d, f in os.walk("test_dvf_out/"):
    for filename in f:
        if filename.endswith("_1kb_all.fasta_gt1000bp_dvfpred.txt"):
            fullfile = os.path.join(r,filename)
            dvf1_files.append(fullfile)
dvf1df = pd.concat([pd.read_csv(fp, sep="\t",low_memory=False) for fp in dvf1_files])
dvf1_cols = []
for y in dvf1df.columns:
	if y == "name":
		dvf1_cols.append("contig")
	else:
		newy = "dvf_" + y
		dvf1_cols.append(newy)

dvf1df.columns = dvf1_cols



merge1 = pd.merge(sag_lens,vs2df,on="contig",how="left")
merge1.fillna(0,inplace=True)
merge2 = pd.merge(merge1,dvf1df,on="contig",how="left")
merge2.fillna(0,inplace=True)



checkv_files = []
trimseq2len = {}
for r, d, f in os.walk("test_checkv_out/"):
    for filename in f:
        if filename.startswith("quality_summary.t"):
            fullfile = os.path.join(r,filename)
            checkv_files.append(fullfile)
        elif filename.endswith("viruses.fna"):
            seqfile = os.path.join(r,filename)
            for record in SeqIO.parse(seqfile,'fasta'):
                name = record.id
                length = len(record.seq)
                trimseq2len[name] = length





checkv_df = pd.concat([pd.read_csv(fp, sep="\t",low_memory=False) for fp in checkv_files])
checkv_cols = checkv_df.columns
new_checkv_cols = []
for p in checkv_cols:
    if p == "contig_id":
        new_checkv_cols.append("contig")
    else:
        newp = p + "_checkv"
        new_checkv_cols.append(newp)
checkv_df.columns = new_checkv_cols






merge3 = pd.merge(merge2,checkv_df,on="contig",how="left")
merge3.fillna(0,inplace=True)


checkvlens = pd.DataFrame.from_dict(trimseq2len,orient='index')
checkvlens.reset_index(inplace=True)
checkvlens.columns = ["contig","checkv_trim_len"]
merge4 = pd.merge(merge3,checkvlens,on="contig",how="left")
merge4.fillna(0,inplace=True)


def get_viral_status(indf):
    virstat = str()
    if indf["vs2_max_score"] >= 0.9:
        virstat = "viral"
    elif indf["vs2_max_score"] >= 0.5 and indf["vs2_hallmark"] >= 2:
        virstat = "viral"
    elif indf["dvf_score"] >= 0.5 and indf["dvf_pvalue"] < 0.05:
        virstat = "viral"
    else:
        virstat = "nonviral"
    return virstat


merge4["viral_status"] = merge4.apply(get_viral_status,axis=1)

merge4.to_csv("test_contig_viral_results_merge.tsv",sep="\t",index=False)

vircon_count = merge4.groupby(["genome","viral_status"]).size()
vircon_count = vircon_count.to_frame()
vircon_count.reset_index(inplace=True)
vircon_count.columns = ["genome","viral_status","count"]
vircon_wide = vircon_count.pivot(index='genome',columns='viral_status',values='count')
vircon_wide.fillna(0,inplace=True)
vircon_wide.reset_index(inplace=True)
vircon_wide.columns = ["genome","nonviral_contigs","viral_contigs"]

vircon_len = merge4.groupby(["genome","viral_status"])["length"].sum()
vircon_len = vircon_len.to_frame()
vircon_len.reset_index(inplace=True)
vircon_len.columns = ["genome","viral_status","sumlen"]
virlen_wide = vircon_len.pivot(index='genome',columns='viral_status',values='sumlen')
virlen_wide.fillna(0,inplace=True)
virlen_wide.reset_index(inplace=True)
virlen_wide.columns = ["genome","viral_length","nonviral_length"]


sagmerge1 = vircon_wide.merge(virlen_wide,on=["genome"],how="left")
sagmerge1.to_csv("test_sag_virinfo.tsv",sep="\t",index=False)



