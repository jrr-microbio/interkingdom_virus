#! /usr/bin/env python
import os
import sys
import pysam
import pybedtools
import argparse
import csv
import pandas as pd
import numpy as np
import math


def fetch_arguments(parser):
	parser.add_argument('--bam_file','-b', dest='bam_file', required=True, default='none',help='input bam file')
	parser.add_argument('--gff_file','-g', dest='gff_file', required=True, default='none',help='input gff file (can also be tab output of geNomad)')
	parser.add_argument('--out_dir','-o', dest='out_dir', required=False, default='./',help='output directory (default = ./)')
	parser.add_argument('--tmp_dir','-t', dest='tmp_dir', required=False, default='tmp/',help='temp file directory (default = tmp/')
	parser.add_argument('--num_threads','-p', dest='threads', required=False, default=1,help='number of threads (default = 1)')
	parser.add_argument('--force_redo', dest='force', required=False, default=False, action="store_true",help='to force full recompute, even if files already exist')


def get_genome_list(args):
	genome_file = args['genome_file']
	print(f"### Creating genome file {genome_file} to help bedtools")
	bam_file = pysam.AlignmentFile(args['bam_file'], "rb")
	header = bam_file.header.as_dict()["SQ"]
	with open(genome_file, 'w') as out_file:
		outwriter = csv.writer(out_file, delimiter='\t',quotechar='"', quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
		for key in header:
			outwriter.writerow([key['SN'],key['LN']])

def tab_to_bed(args):
	df_genes = pd.read_csv(args['gff_file'], sep='\t', low_memory=False)
	df_genes["contig"] = df_genes["gene"].apply(lambda x: ("_").join(x.split("_")[:-1])) ## Guess the contig name
	df_genes["strand"] = df_genes["strand"].apply(lambda x: clean_strand(str(x))) ## Clean up the strand
	## Bed file number sequence from 0 to n (i.e. first base is coordinate 0)
	## When defining feature in a bed file, the end coordinate is not included, i.e. a bed file with columns #2 and #3 or 0 and 100 would only reference to positions 0 to 99
	## Since IMG and other pipelines defines coordinates starting from 1, that means we need to decrease start by 1 to switch to 0-based coordinates AND include an extra base in 3'
	## For more information, see https://genome.ucsc.edu/FAQ/FAQformat.html
	df_genes["start"] = df_genes["start"] - 1 ## Need to decrease start by 1, see above
	df_genomes = pd.read_csv(args['genome_file'], header=None, names=['genome','length'], sep='\t')
	df_genes["contig"] = pd.Categorical(df_genes['contig'], df_genomes['genome'], ordered=True) ## prepare to sort according to the bam file
	if df_genes["contig"].isnull().any():
		print("STOPPED -- SOME SEQUENCES WERE IN THE ANNOTATION FILE BUT NOT IN THE BAM FILE")
		## Should we handle these cases ? will probably fail at bedtools step anyway
		sys.exit(1)
	df_genes = df_genes.sort_values(by = ["contig","start","gene"], ascending = [True,True,True])
	df_genes['bed_score'] = 0
	df_genes['bed_type'] = 'CDS'
	df_genes = df_genes.rename(columns = {'annotation_description' : 'product'} )
	df_genes[['contig','start','end','gene','bed_score','strand','bed_type','product']].to_csv(args['bed_file'], sep='\t', header=None, index=False, na_rep='NA')


def clean_strand(st):
	##
	if st == "1":
		st = "+"
	if st == "-1":
		st = "-"
	# if st == "+":
	# 	st = "1"
	# if st == "-":
	# 	st = "-1"
	return(st)

def gff_to_bed(args):
	store = {}
	print(f"Reading {args['gff_file']}")
	with open(args['gff_file'], "r") as f:
		reader = csv.reader(f, delimiter="\t")
		i = 0
		for row in reader:
			if not(row[0].startswith("#")):
				# print(row)
				## We assume standard gff
				contig = str(row[0])
				bed_type = str(row[2])
				## Bed file number sequence from 0 to n (i.e. first base is coordinate 0)
				## When defining feature in a bed file, the end coordinate is not included, i.e. a bed file with columns #2 and #3 or 0 and 100 would only reference to positions 0 to 99
				## Since IMG and other pipelines defines coordinates starting from 1, that means we need to decrease start by 1 to switch to 0-based coordinates AND include an extra base in 3'
				## For more information, see https://genome.ucsc.edu/FAQ/FAQformat.html
				start = int(row[3]) - 1 ## Need to decrease start by 1, see above
				end = int(row[4])
				strand = clean_strand(str(row[6])) ## We convert the strand from 1 and -1 instead tp + and - (official bed format)
				tab_h = row[8].split(";")
				gene_id = "NA"
				product = "NA"
				bed_score = 0
				for cell in tab_h:
					if cell.startswith("ID="):
						gene_id=str(cell[3:])
					if cell.startswith("product="):
						product=str(cell[8:])
				if contig not in store:
					store[contig] = {}
				# print(f"{gene_id} // {product}")
				i = i+1 ## We assume the gff is sorted by start position
				if bed_type!="CRISPR" and bed_type!="repeat_unit": ## Ignore problematic CRISPR / CRISPR arrays
					store[contig][i] = [contig,start,end,gene_id,bed_score,strand,bed_type,product]
					# print(f"Loaded one more gene for contig {contig} ({i}) // {store[contig][i]}")
					# 'contig','start','end','gene','bed_score','strand','bed_type','product'
	print(f"Writing {args['bed_file']} based on the genome order in {args['genome_file']}")
	with open(args['genome_file'], "r") as f, open(args['bed_file'], 'w') as o:
		reader = csv.reader(f, delimiter="\t")
		writer = csv.writer(o, delimiter="\t")
		for row in reader:
			# print(f"{row[0]}")
			if row[0] in store:
				for i in store[row[0]]:
					# print(i)
					# print(f"{store[row[0]][i]}")
					writer.writerow(store[row[0]][i])



def split_bam(args):
	print(f"### Split the bam file {args['bam_file']} into {args['final_f']} and {args['final_r']}")
	## Note - some of the magic used with pysam is from https://github.com/pysam-developers/pysam/issues/677
	sorted_bam = os.path.join(args['tmp_dir'],"entire_sorted.bam")
	## sort and index
	pysam.sort("-@", args['threads'], "-o", sorted_bam, args['bam_file'])
	pysam.index("-@", args['threads'], sorted_bam)
	## get only the forward reads, merge and index
	tmp_bam_f1 = os.path.join(args['tmp_dir'],"tmpfwd1.bam")
	fh = open(tmp_bam_f1, 'w')
	fh.close()
	tmp_bam_f2 = os.path.join(args['tmp_dir'],"tmpfwd2.bam")
	fh = open(tmp_bam_f2, 'w')
	fh.close()
	## NOTE: Explanation of bitwise flag available at https://broadinstitute.github.io/picard/explain-flags.html
	# Forward strand.
	# 1. alignments of the second in pair if they map to the forward strand
	# 2. alignments of the first in pair if they map to the reverse  strand
	pysam.view("-@", args['threads'], "-b", "-f", "128", "-F", "16", sorted_bam, "-o", tmp_bam_f1, save_stdout=tmp_bam_f1, catch_stdout=False)
	print(f"{tmp_bam_f1} created")
	pysam.view("-@", args['threads'], "-b", "-f", "80", sorted_bam, "-o", tmp_bam_f2, save_stdout=tmp_bam_f2, catch_stdout=False)
	print(f"{tmp_bam_f2} created")
	# Combine alignments that originate on the forward strand.
	pysam.merge("-@", args['threads'], "-f", args['final_f'], tmp_bam_f1, tmp_bam_f2)
	pysam.index("-@", args['threads'], args['final_f'])
	print(f"{args['final_f']} created and indexed")
	#
	# Reverse strand
	tmp_bam_r1 = os.path.join(args['tmp_dir'],"tmprev1.bam")
	fh = open(tmp_bam_r1, 'w')
	fh.close()
	tmp_bam_r2 = os.path.join(args['tmp_dir'],"tmprev2.bam")
	fh = open(tmp_bam_r2, 'w')
	fh.close()
	# 1. alignments of the second in pair if they map to the reverse strand
	# 2. alignments of the first in pair if they map to the forward strand
	pysam.view("-@", args['threads'], "-b", "-f", "144", sorted_bam, "-o", tmp_bam_r1, save_stdout=tmp_bam_r1, catch_stdout=False)
	print(f"{tmp_bam_r1} created")
	pysam.view("-@", args['threads'], "-b", "-f", "64", "-F", "16", sorted_bam, "-o", tmp_bam_r2, save_stdout=tmp_bam_r2, catch_stdout=False)
	print(f"{tmp_bam_r2} created")
	# Combine alignments that originate on the reverse strand.
	pysam.merge("-@", args['threads'], "-f", args['final_r'], tmp_bam_r1, tmp_bam_r2, catch_stdout=False)
	pysam.index("-@", args['threads'], args['final_r'])
	print(f"{args['final_r']} created and indexed")
	if os.path.exists(args['final_f']) and os.path.exists(args['final_r']):
		os.remove(tmp_bam_f1)
		os.remove(tmp_bam_f2)
		os.remove(tmp_bam_r1)
		os.remove(tmp_bam_r2)
	else:
		print("we had some issue when trying to split the bam file")
		sys.exit(1)

def stat_cover(x):
	# print(f"{x}")
	# print(f"{x['depth']}")
	# print(f"{x['bases']}")
	avg = np.average(x['depth'], weights = x['bases'])
	variance = np.average((x['depth']-avg)**2, weights=x['bases'])
	stdev = math.sqrt(variance)
	# median = np.quantile(x['depth'], 0.5, method='inverted_cdf', weights=x['bases']) ## Need numpy 2.0.0, not yet fully in conda
	median = weighted_quantiles(np.array(x['depth']), np.array(x['bases']), 0.5, False)
	# if avg > 0:
		# print(f"{x}")
		# print(f"{{avg} // {stdev} // {median}")
	return pd.Series([avg, stdev, median],index=['average','stdev','median'])

## Taken from https://stackoverflow.com/questions/20601872/numpy-or-scipy-to-calculate-weighted-median
def weighted_quantiles(values, weights, quantiles=0.5, interpolate=False):
	i = values.argsort()
	sorted_weights = weights[i]
	sorted_values = values[i]
	Sn = sorted_weights.cumsum()
	if interpolate:
		Pn = (Sn - sorted_weights/2 ) / Sn[-1]
		return np.interp(quantiles, Pn, sorted_values)
	else:
		# print("We do not interpolate")
		# print(f"We will search {Sn} for {quantiles}")
		return sorted_values[np.searchsorted(Sn, quantiles * Sn[-1])]

def get_coverage_stats(bam_file,b,args):
	print(f"## Reading {bam_file}")
	x = pybedtools.BedTool(bam_file)
	c = b.coverage(x,sorted=True,g=args['genome_file'])
	cover = {}
	print("first reading total reads and fraction covered per gene")
	df = c.to_dataframe()
	# chrom start end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts
	## The above are the default column names for a bed file, we replace when needed
	df.rename(columns={"chrom": "scaffold_id", "name": "gene_id"}, inplace=True)
	if "thickStart" in df.columns:
		df.rename(columns={"thickStart": "locus_type"}, inplace=True)
	if "thickEnd" in df.columns:
		df.rename(columns={"thickEnd": "product"}, inplace=True)
	## The last four columns are always reads / base / size / fraction covered, so we rename as well
	## Specifically, bedtools adds four columns at the end (we use indexing from the end to get these last four columns):
	# -4 - The number of features (reads) in B that overlapped (by at least one base pair) the A interval.
	# -3 - The number of bases in A that had non-zero coverage from features (reads) in B.
	# -2 - The length of the entry in A.
	# -1 - The fraction of bases in A that had non-zero coverage from features in B.
	df.rename(columns={df.columns[-4]: "total_reads", df.columns[-3]: "total_bases", df.columns[-2]: "length", df.columns[-1]: "fraction"}, inplace=True)
	# print(f"{df}")
	print("now getting the average, median, and stdev of read depth")
	chist = b.coverage(x,sorted=True,g=args['genome_file'],hist=True)
	## This time the last four columns are still depth / base / size / fraction but because we asked for "hist=True", we'll get the actual distribution and counts
	## Specifically, bedtools adds four columns at the end (we use indexing from the end to get these last four columns):
	# -4 - The number of features in B that overlapped (by at least one base pair)
	# -3 - The number of bases with the corresponding number of features in B that are overlapped
	# -2 - The length of the entry in A.
	# -1 - The fraction of bases in A that had the corresponding overlap from features in B.
	dfhist = chist.to_dataframe()
	dfhist = dfhist[dfhist['chrom']!="all"]
	dfhist.rename(columns={"chrom": "scaffold_id", "name": "gene_id"}, inplace=True)
	dfhist.rename(columns={dfhist.columns[-4]: "depth", dfhist.columns[-3]: "bases", dfhist.columns[-2]: "length", dfhist.columns[-1]: "fraction"}, inplace=True)
	# print(f"{dfhist}")
	dfhist_stats = dfhist.groupby(['scaffold_id','gene_id']).apply(lambda x: stat_cover(x), include_groups=False).reset_index()
	print(f"{dfhist_stats}")
	## Now merging with the first df
	df_final = df.merge(dfhist_stats,left_on=['scaffold_id','gene_id'],right_on=['scaffold_id','gene_id'],how="left")
	print(f"{df_final}")
	return(df_final)

def write_outfile(args,df_cover):
	print(f"{df_cover}")
	df_cover.to_csv(args['out_file'], index=False, sep='\t')
	## And do some stat as QC
	tr_f = df_cover.groupby(['strand'])['reads_cnt_plus'].apply(lambda x: sum(x))
	print(f"Total reads on forward {tr_f}")
	tr_r = df_cover.groupby(['strand'])['reads_cnt_minus'].apply(lambda x: sum(x))
	print(f"Total reads on reverse {tr_r}")
	tr = tr_f["+"] + tr_f["-"] + tr_r["+"] + tr_r["-"]
	tr_expected = tr_f["+"] + tr_r["-"]
	if tr>0:
		p_exp = tr_expected / tr * 100
	else:
		p_exp = NA
	print(f"That gives {tr} total reads and {tr_expected} mapping on the expected strand, i.e. ({p_exp:0.2f} %) expected mapping rate")
	## Checking that we recovered all the genes too
	df_bed = pd.read_csv(args['bed_file'], sep='\t', low_memory=False, header=None)
	tg_bed = len(df_bed)
	print(f"{df_bed}")
	tg_cover = len(df_cover)
	print(f"{df_cover}")
	with open(args['out_stat_log'], 'w') as stat_file:
		stat_file.write(f"Total genes in bed file\t{tg_bed}\n")
		stat_file.write(f"Total genes in count file\t{tg_cover}\n")
		if tg_bed != tg_cover:
			print(f"#!#!#!#!#!#! WARNING #!#!#!#!#!#!#")
			print(f"We did not find the same number of genes between the bed file {tg_bed} and the count file {tg_cover}, something may be wrong")
			stat_file.write("WARNING - INCONSISTENT NUMBER OF GENES\n")
		stat_file.write(f"Total reads mapped to genes\t{tr}\n")
		stat_file.write(f"Expected reads for genes on plus strand\t{tr_f['+']:,}\n")
		stat_file.write(f"Unexpected reads for genes on plus strand\t{tr_f['-']:,}\n")
		stat_file.write(f"Expected reads for genes on minus strand\t{tr_r['-']:,}\n")
		stat_file.write(f"Unexpected reads for genes on minus strand\t{tr_r['+']:,}\n")
		stat_file.write(f"Overall percentage of expected mapping\t{p_exp:0.2f} %\n")

def main():
	parser = argparse.ArgumentParser()
	fetch_arguments(parser)
	args = vars(parser.parse_args())
	## First verify that the input bam file is here
	if not os.path.exists(args["bam_file"]):
		print(f"Pblm, I could not find file {args['bam_file']}, and I do need an input file")
		sys.exit(1)
	if not os.path.isdir(args["out_dir"]):
		print(f"Pblm, I could not find file {args['out_dir']}, and I do need an output folder")
		sys.exit(1)
	if not os.path.isdir(args["tmp_dir"]):
		print(f"Pblm, I could not find file {args['tmp_dir']}, and I do need an output folder")
		sys.exit(1)
	### Create an ordered list of genome from original bam file
	args['genome_file'] = os.path.join(args['tmp_dir'],"genomes.txt")
	if (os.path.exists(args['genome_file']) and not(args['force'])):
		print(f"{args['genome_file']} already here, we skip")
	else:
		get_genome_list(args)
	### Read gff (or geNomad tab) and create a sorted bed file instead (sorted based on the order of the contigs in the bam file)
	args['bed_file'] = os.path.join(args['tmp_dir'],'genes.bed')
	if (os.path.exists(args['bed_file']) and not(args['force'])):
		print(f"{args['bed_file']} already here, we skip")
	else:
		if args['gff_file'].endswith("_genes.tsv"):
			print(f"We guess the file {args['gff_file']} is actually the output tab file from geNomad, so we will use the tab_to_bed conversion")
			tab_to_bed(args)
		else:
			print(f"We guess the file {args['gff_file']} is a standard gff file, so we will use the gff_to_bed conversion")
			gff_to_bed(args)
	### Use samtools to obtain two files:
	# First one with all cases in which the original fragment was from the forward strand (F2 - R1), i.e. take all reads which are first in pair and on reverse, and all reads which are second in pair and on forward
	# Second one with the opposite, i.e. all cases in which the original fragment was from the reverse strand (F1 - R2)
	args['final_f'] = os.path.join(args['tmp_dir'],'fwd_sorted.bam');
	args['final_r'] = os.path.join(args['tmp_dir'],'rev_sorted.bam');
	if (os.path.exists(args['final_f']) and os.path.exists(args['final_r']) and not(args['force'])):
		print(f"{args['final_f']} and {args['final_r']} already here, we skip")
	else:
		split_bam(args)
	### Use bedtools to get a distribution of read counts for each feature, and calculate from this the average, median, and standard deviation
	# Load bed file information
	b = pybedtools.BedTool(args['bed_file'])
	# Get the stats
	print("### Get coverage stats")
	print("first reads mapped on forward strand")
	cover_f = get_coverage_stats(args['final_f'],b,args);
	## Rename columns in cover_f and select only the relevant columns
	cover_f.rename(columns={cover_f.columns[-7]: "reads_cnt_plus", cover_f.columns[-4]: "covered_fraction_plus", cover_f.columns[-3]: "mean_plus", cover_f.columns[-2]: "stdev_plus", cover_f.columns[-1]: "median_plus"}, inplace=True)
	cover_f = cover_f[['gene_id', 'scaffold_id', 'strand', 'length', 'locus_type', 'reads_cnt_plus', 'mean_plus', 'median_plus', 'stdev_plus', 'covered_fraction_plus']]
	print("then reads mapped on reverse strand")
	cover_r = get_coverage_stats(args['final_r'],b,args);
	## Rename columns in cover_r and select only the relevant columns
	cover_r.rename(columns={cover_r.columns[-7]: "reads_cnt_minus", cover_r.columns[-4]: "covered_fraction_minus", cover_r.columns[-3]: "mean_minus", cover_r.columns[-2]: "stdev_minus", cover_r.columns[-1]: "median_minus"}, inplace=True)
	if "product" in cover_r:
		cover_r = cover_r[['gene_id', 'scaffold_id', 'reads_cnt_minus', 'mean_minus', 'median_minus', 'stdev_minus', 'covered_fraction_minus', 'product']]
	else:
		cover_r = cover_r[['gene_id', 'scaffold_id', 'reads_cnt_minus', 'mean_minus', 'median_minus', 'stdev_minus', 'covered_fraction_minus']]
	cover_merged = cover_f.merge(cover_r, left_on=['gene_id','scaffold_id'], right_on=['gene_id','scaffold_id'], how='outer', suffixes=(False, False))
	### Prepare the output file
	args["out_file"] = os.path.join(args["out_dir"],"rnaseq_gea.tab")
	args["out_stat_log"] = os.path.join(args["out_dir"],"rnaseq_stats.log")
	print(f"### Write the final output file {args['out_file']} as well as the summary stats file {args['out_stat_log']}")
	write_outfile(args,cover_merged)
	print(f"###### All done ! ######")





if __name__ == "__main__":
	output = main()
