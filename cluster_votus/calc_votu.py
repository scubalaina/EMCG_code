import os, sys, re, argparse
import pandas as pd

args_parser = argparse.ArgumentParser(description="Script clustering genomes into vOTUs based on ANI", epilog="Bigelow Laboratory for Ocean Sciences")
args_parser.add_argument('-i', '--input', required=True, help='BLASTn output of genomes against each other in output format 6.')
args_parser.add_argument('-l', '--lens', required=True, help='Tab-seperated file with each open reading frame length of all genomes in BLAST search. First column is gene; second column is length.')
args_parser.add_argument('-n', '--name', required=True, help='Label for output merged files with suffix _votu_membership.tsv')

args_parser = args_parser.parse_args()

blast = args_parser.input
orffile = args_parser.lens
name = args_parser.name

votu_mem = name + "_votu_membership.tsv"


blastdf = pd.read_csv(blast,sep="\t",header=None,index_col=False)

blastdf.columns = ["geneA","geneB","id","alnlen","mismatch","gapopen","qstart","qend", "sstart","send","bitscore","evalue"]
# take the best alignment of two genes to each other based on bitscore 
bestdf = blastdf.sort_values(["bitscore"],ascending=False).drop_duplicates(subset =["geneA","geneB"], keep="first")
# make genome A column (genA) and genome B column (genB)
bestdf["genA"] = bestdf["geneA"].str.split("_").str[0:2].str.join("_")
bestdf["genB"] = bestdf["geneB"].str.split("_").str[0:2].str.join("_")

# remove self hits
#nonself = bestdf[bestdf["genA"] != bestdf["genB"]]

# get a gene's best hit to another genome based on bitscore
best_genhit = bestdf.sort_values(["bitscore"],ascending=False).drop_duplicates(subset=["geneA","genB"],keep="first")

# count number of unique hits from genA to genB
uniqhits = best_genhit.groupby(["genA","genB"]).size()
uniqhits = uniqhits.to_frame()
uniqhits.reset_index(inplace=True)
uniqhits.columns = ["genA","genB","genA_hits"]

# get ANI 
ani = best_genhit.groupby(["genA","genB"])["id"].mean()
ani = ani.to_frame()
ani.reset_index(inplace=True)
ani.columns = ["genA","genB","genA_ani"]

# get aligned length of genA to genB
best_genhit["genA_alnlen"] = best_genhit["qend"] - best_genhit["qstart"]
best_genhit["genA_alnlen"] = best_genhit["genA_alnlen"].abs()

alnlen = best_genhit.groupby(["genA","genB"])["genA_alnlen"].sum()
alnlen = alnlen.to_frame()
alnlen.reset_index(inplace=True)
alnlen.columns = ["genA","genB","genA_sumaln"]

# merge parsed outputs to get initial ani and coverage table
merge1 = uniqhits.merge(ani,on=["genA","genB"],how="left")
merge2 = merge1.merge(alnlen,on=["genA","genB"],how="left")

orflens = pd.read_csv(orffile,sep="\t",index_col=False,header=0)
orflens["genome"] = orflens["gene"].str.split("_").str[0:2].str.join("_")

genlens = orflens.groupby(["genome"]).size()
genlens = genlens.to_frame()
genlens.reset_index(inplace=True)
genlens.columns = ["genA","genA_genes"]

sumlen = orflens.groupby(["genome"])["length"].sum()
sumlen = sumlen.to_frame()
sumlen.reset_index(inplace=True)
sumlen.columns = ["genA","genA_cdslen"]

merge3 = merge2.merge(genlens,on=["genA"],how="left")
merge4 = merge3.merge(sumlen,on=["genA"],how="left")

merge4["prop_genA_cdslen"] = merge4["genA_sumaln"] / merge4["genA_cdslen"]
merge4["prop_genA_genes"] = merge4["genA_hits"] / merge4["genA_genes"]
votu1 = merge4[merge4["prop_genA_cdslen"] >= 0.85]
votu2 = votu1[votu1["genA_ani"] >= 95]

import networkx as nx
G = nx.from_pandas_edgelist(votu2, source='genA', target='genB')
subgraphs = list(nx.connected_components(G))
gen2pop = {}
poplist = []
for idx, subgraph_nodes in enumerate(subgraphs, start=1):
    subgraph = G.subgraph(subgraph_nodes)
    popname = "vOTU" + "_" + str(idx)
    for x in subgraph_nodes:
        #print(popname + "\t" + x)
        gen2pop[x] = popname
    poplist.append(popname)

popdf = pd.DataFrame.from_dict(gen2pop,orient='index')
popdf.reset_index(inplace=True)
popdf.columns = ['genome',"vOTU"]


#popdf.to_csv(votu_mem,sep="\t",index=False)


