import os, sys, re, argparse
import pandas as pd
from collections import defaultdict

args_parser = argparse.ArgumentParser(description="Script for running prodigal on SAGs", epilog="Bigelow Laboratory for Ocean Sciences")
args_parser.add_argument('-i', '--infile', required=True, help='Input directory of parsed HMM tables with suffix _nog_parsed.tsv .')
args_parser.add_argument('-l', '--label', required=True, help='Label for output merged files with suffixes _nog_parsed.tsv and _nog_parsed_minbit30.tsv.')

args_parser = args_parser.parse_args()

infile = args_parser.infile
label = args_parser.label

indf = pd.read_csv(infile,sep="\t",header=0,index_col=False)
outfile = label + "_domain_count.tsv"

sagdomcount = indf.groupby(["genome","domain"]).size()
sagdomcount = sagdomcount.to_frame()
sagdomcount.reset_index(inplace=True)
sagdomcount.columns = ["genome","domain","genes"]
sagdom_wide = sagdomcount.pivot(values = "genes", index=["genome"], columns = ["domain"])
sagdom_wide.reset_index(inplace=True)
sagdom_wide["cellular_hits"]  = sagdom_wide["Archaea"] + sagdom_wide["Bacteria"] + sagdom_wide["Eukarya"]
sagdom_wide.to_csv(outfile,sep="\t",index=False)

