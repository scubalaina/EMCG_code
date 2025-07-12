import os, sys, re, argparse
import pandas as pd

args_parser = argparse.ArgumentParser(description="Script for categorizing SAG entities", epilog="Bigelow Laboratory for Ocean Sciences")
args_parser.add_argument('-i', '--info', required=True, help='Input directory of parsed HMM tables with suffix _nog_parsed.tsv .')
args_parser.add_argument('-d', '--domfile', required=True, help='TSV file where each row is a SAG and the columns are the number of genes in each domain.')
args_parser.add_argument('-v', '--virfile', required=True, help='TSV file where each row is a SAG and the columns are the viral content of that SAG')

args_parser = args_parser.parse_args()

info = args_parser.info
domfile = args_parser.domfile
virfile = args_parser.virfile

outfile = re.sub(".tsv","_classified.tsv",info)

infodf = pd.read_csv(info,sep="\t",index_col=False,header=0)
domdf = pd.read_csv(domfile,sep="\t",index_col=False,header=0)
virdf = pd.read_csv(virfile,sep="\t",index_col=False,header=0)

merge1 = infodf.merge(domdf,on=["genome"],how="left")
merge1["prop_cellgen"] = merge1["cell_gen"] / merge1["num_genes"]

merge2 = merge1.merge(virdf,on=["genome"],how="left")
merge2["prop_virlen"] = merge2["viral_length"] / merge2["untrimmed_length"]


def get_entity(indf):
	entity = str()
	if indf["prop_cellgen"] > indf["prop_virlen"]:
		entity = "cellular"
	elif indf["prop_cellgen"] < indf["prop_virlen"]:
		entity = "viral"
	else:
		entity = "undetermined"
	return entity

merge2["entity"] = merge2.apply(get_entity,axis=1)
merge2.to_csv(outfile,sep="\t",index=False)



