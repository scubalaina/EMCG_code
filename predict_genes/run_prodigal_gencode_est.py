import os, re, sys, shlex, subprocess, shutil, argparse
from Bio import SeqIO

args_parser = argparse.ArgumentParser(description="Script for running prodigal on SAGs", epilog="Bigelow Laboratory for Ocean Sciences")
args_parser.add_argument('-i', '--indir', required=True, help='Input directory of individual genome nucleotide fasta files ending with .fna or .fasta')
args_parser.add_argument('-o', '--outdir', required=True, help='Output directory where amino acid sequences (.faa) and gene nucleotide sequences (orf_nucl.fna) are deposited.')
args_parser.add_argument('-f', '--code_file', required=True, help='Tab-seperated file with each genome, its best genetic code, its coding density, and its average gene length.')

args_parser = args_parser.parse_args()

indir = args_parser.indir
outdir = args_parser.outdir
outfile_name = args_parser.code_file

if os.path.isdir(outdir):
	pass
else:
	os.mkdir(outdir)

if indir.endswith("/"):
	pass
else:
	indir = indir + "/"

if outdir.endswith("/"):
	pass
else:
	outdir = outdir + "/"

indir = sys.argv[1]
outdir = sys.argv[2]
outfile = open(outfile_name,'a')
codes = ['1', '2', '3', '4', '5', '6', '9', '10', '11', '12', '13', '14', '15', '16', '21', '22', '23', '24', '25']

outfile.write("genome" + "\t" + "code" + "\t" + "coding_density" + "\t" + "avg_protlen" + "\n")

def get_density(protein,nucl_length):
	prot_length = float(0)
	prot_lens = []
	for j in SeqIO.parse(protein, "fasta"):
		prot_length += len(j.seq)
		orf_len = len(j.seq) * 3
		prot_lens.append(orf_len)
	density = (prot_length * 3) / nucl_length
	avg_len = sum(prot_lens) / len(prot_lens)
	return density, avg_len 

def get_bestcode(code2density,code2avglen,outdir,genome):
	standard = float(code2density['11'])
	max_key = str()
	if standard < float(0.8):
		max_key = max(code2density, key=code2density. get)		
	else:
		max_key = '11'
	max_avglen = code2avglen[max_key]
	new_protein = outdir + genome + "." + max_key + ".faa"
	new_gene = outdir + genome + "." + max_key + ".orf_nucl.fna"
	max_density = code2density[max_key]
	return max_key, max_avglen, max_density, new_protein, new_gene


for i in os.listdir(indir):
	if i.endswith(".fna") or i.endswith(".fasta"):
		fasta = os.path.join(indir, i)
		nucl_length = float(0)
		#genome = re.sub("_contigs.fasta","",i)
		genome = re.sub(".fasta","",i)
		#genome = re.sub(".fna","",i)

		for j in SeqIO.parse(fasta, "fasta"):
			nucl_length += len(j.seq)

		code2density = {}
		code2avglen = {}
		for code in codes:
			#protein = os.path.join(outdir, re.sub(".fna", "."+str(code)+".faa", i))
			protein = outdir + genome + "." + code + ".faa"
			genes = outdir + genome + "." + code + ".orf_nucl.fna"
			gff = outdir + genome + "." + code + ".gff"
			#print(protein)

			if nucl_length < 20000:
				cmd = "prodigal -p meta -i "+ fasta +" -a "+ protein + " -d " + genes + " -g "+ code + " -o " + gff

			else:
				cmd = "prodigal -i "+ fasta +" -a "+ protein + " -d " + genes +" -g "+ code + " -o " + gff

			cmd2 = shlex.split(cmd)
			#print(cmd)
			subprocess.call(cmd2, stdout=open("out.txt", "w"), stderr=open("err.txt", "w"))

			
			file_size = os.path.getsize(protein)
			if file_size > 0:
				density, avglen = get_density(protein,nucl_length)
				code2density[code] = density
				code2avglen[code] = avglen
		
		if len(code2density) > 0:
			#print(genome)
			max_key, max_avglen, max_density, new_protein, new_gene = get_bestcode(code2density,code2avglen,outdir,genome)
			bestcode_list = [genome,max_key,str(max_density),str(max_avglen)]
			bestcode_info = "\t".join(bestcode_list)
			#print(bestcode_info)
			outfile.write(bestcode_info + "\n")
			for filenameX in os.listdir(outdir):
				if filenameX.startswith(genome):
					fullfileX = outdir + filenameX
					if fullfileX == new_gene:
						pass
					elif fullfileX == new_protein:
						pass
					else:
						os.remove(fullfileX)
			#outfile.write(bestcode_info + "\n")
		

for filename in os.listdir(outdir):
	if filename.endswith(".faa"):
		fullfile = outdir + filename
		filelist = filename.split(".")
		prefixlist = filelist[0:-2]
		prefix = ".".join(prefixlist)
		suffix = filelist[-1]
		newname = prefix + "." + suffix
		newout = outdir + newname
		shutil.move(fullfile,newout)
	elif filename.endswith("orf_nucl.fna"):
		fullfile = outdir + filename
		filelist = filename.split(".")
		prefixlist = filelist[0:-3]
		suffixlist = filelist[-2:]
		prefix = ".".join(prefixlist)
		suffix = ".".join(suffixlist)
		newname = prefix + "_" + suffix
		newout = outdir + newname
		shutil.move(fullfile,newout)
	elif filename.endswith(".gff"):
		fullfile = outdir + filename
		filelist = filename.split(".")
		prefixlist = filelist[0:-2]
		suffixlist = filelist[-1:]
		prefix = ".".join(prefixlist)
		suffix = ".".join(suffixlist)
		newname = prefix + "_" + suffix
		newout = outdir + newname
		shutil.move(fullfile,newout)



