import os, sys, re, argparse
import pandas as pd

args_parser = argparse.ArgumentParser(description="Script for running prodigal on SAGs", epilog="Bigelow Laboratory for Ocean Sciences")
args_parser.add_argument('-i', '--indir', required=True, help='Input directory of parsed HMM tables with suffix _nog_parsed.tsv .')
args_parser.add_argument('-l', '--label', required=True, help='Label for output merged files with suffixes _nog_parsed.tsv and _nog_parsed_minbit30.tsv.')

args_parser = args_parser.parse_args()

indir = args_parser.indir
dataset = args_parser.label


parsed = []
for r,d,f in os.walk(indir):
	for filename in f:
		if filename.endswith("nog_parsed.tsv"):
			fullfile = os.path.join(r,filename)
			parsed.append(fullfile)


df = pd.concat([pd.read_csv(fp, sep="\t",low_memory=False)for fp in parsed])
df['genome'] = df['protein'].str.split('_').str[0:2].str.join("_")
sys_first_column = df.pop('genome')
df.insert(0, 'genome', sys_first_column)
outfilename = dataset + "_nog_parsed.tsv"
overview_open = open(outfilename,'w')
df.to_csv(overview_open,sep="\t",index=False)

df_bit30 = df[df['bitscore'] >= 30]
out30 = re.sub(".tsv","_minbit30.tsv",outfilename)
df_bit30.to_csv(out30,sep="\t",index=False)
