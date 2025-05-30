import os, sys, re, argparse
from collections import defaultdict

#  python parse_nog_hmm.py -i microcap_capsuleA1_nog_out/ -o microcap_capsuleA1_nog_parsedvir/ -n nog_annotation_virupdated.tsv

args_parser = argparse.ArgumentParser(description="Script for running prodigal on SAGs", epilog="Bigelow Laboratory for Ocean Sciences")
args_parser.add_argument('-i', '--indir', required=True, help='Input directory of individual tab-seperated HMMsearch output tables file with _hmm_tbl.txt suffix')
args_parser.add_argument('-o', '--outdir', required=True, help='Output directory where parsed HMM tables are deposited with _nog_parsed.tsv suffix.')
args_parser.add_argument('-n', '--nog_meta', required=True, help='Metadata file with NOG HMM profile information called nog_annotation_virupdated.tsv')

args_parser = args_parser.parse_args()

infolder = args_parser.indir
outfolder = args_parser.outdir
nog_meta = args_parser.nog_meta

annot = open(nog_meta,'r')

acc2cat = {}
acc2name = {}
acc2domain = {}
acc2virdom = {}

EVALUE = 0.00001



if os.path.isdir(outfolder):
	pass
else:
	os.mkdir(outfolder)

if infolder.endswith("/"):
	pass
else:
	infolder = infolder + "/"

if outfolder.endswith("/"):
	pass
else:
	outfolder = outfolder + "/"


for i in annot:
	line = i.rstrip()
	tabs = line.split("\t")
	acc = tabs[0]
	name = tabs[1]
	cat = tabs[2]
	domain = tabs[3]
	virdom = tabs[4]
	acc2name[acc] = name
	acc2domain[acc] = domain
	acc2cat[acc] = cat
	acc2virdom[acc] = virdom

for filename in os.listdir(infolder):
	genome = re.sub("_hmm_tbl.txt","",filename)
	fullfile = infolder + filename
	print(genome)
	outfile = outfolder + genome + "_nog_parsed.tsv"
	input1 = open(fullfile,'r')
	protein2hit_dict = {}
	protein2bit_dict = {}
	protein2evalue_dict = {}
	protein2name_dict = {}
	protein2cat_dict = {}
	protein2domain_dict = {}
	for i in input1.readlines():
		line = i.rstrip()
		if line.startswith("#"):
			pass
		else:
			newline = re.sub("\s+", "\t", line)
			tabs = newline.split("\t")
			query = tabs[0]
			bit_score = float(tabs[5])
			evalue = float(tabs[4])
			hit_list = tabs[2].split(".")
			hit = hit_list[0]
			if query in protein2bit_dict: # If query is in prtein2bit_dict, it means we have seen this protein before, so now we want to see what bit score it had to the previous hit.
				if protein2bit_dict[query] > float(bit_score):
					pass
				elif evalue > EVALUE:
					pass
				else:
					protein2bit_dict[query] = float(bit_score)
					protein2hit_dict[query] = hit
					protein2evalue_dict[query] = float(evalue)
			elif evalue > EVALUE:
				pass
			else:
				protein2bit_dict[query] = float(bit_score)
				protein2hit_dict[query] = hit
				protein2evalue_dict[query] = float(evalue)
	outopen = open(outfile,'w')
	headers = ["protein","accession","bitscore","evalue","description","category","domain","virdom"]
	headerlist = "\t".join(headers)
	outopen.write(headerlist + "\n")
	for protein, hit in protein2hit_dict.items():
		name = acc2name[hit]
		cat = acc2cat[hit]
		domain = acc2domain[hit]
		virdom = acc2virdom[hit]
		evalue = protein2evalue_dict[protein]
		bit = protein2bit_dict[protein]
		info = [protein,hit,str(bit),str(evalue),name,cat,domain,virdom]
		infolist = "\t".join(info)
		outopen.write(infolist + "\n")


