#!/usr/bin/env python
#
# Read VCF file, assign Typhi genotype codes.
#
# Note the SNP coordinates used here are relative to Typhi CT18 (AL513382) 
#   so the input MUST be a VCF obtained via mapping to this reference sequence.
#
# You MUST tell the script the name of the Typhi CT18 chromosome reference in your VCF 
#   file, this is the entry in the first column (#CHROM) of the data part of the file.
#   This is necessary in case you have mapped to multiple replicons (e.g. chromosome and
#   plasmid) and all the results appear in the VCF file.
#
# E.g. if your VCF file looks like this:
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	2332STDY5573284
#AL513382	744	.	A	G	0.213	.	DP=49;VDB=5.579080e-02;RPB=9.166667e-01;AF1=0.495;AC1=1;DP4=6,24,9,0;MQ=60;FQ=-10.4;PV4=2.4e-05,7.9e-08,1,1	GT:PL:DP:SP:GQ	0/1:17,0,255:39:46:5
#AL513382	762	.	A	G	29	.	DP=43;VDB=1.182296e-01;RPB=1.876799e+00;AF1=0.5;AC1=1;DP4=6,23,12,0;MQ=60;FQ=32;PV4=2.4e-06,1.4e-12,1,1	GT:PL:DP:SP:GQ	0/1:59,0,255:41:56:62
#AL513382	3727	.	A	G	30	.	DP=51;VDB=1.213276e-01;RPB=1.690455e+00;AF1=0.5;AC1=1;DP4=7,26,15,0;MQ=60;FQ=33;PV4=1.6e-07,2.4e-17,1,1	GT:PL:DP:SP:GQ	0/1:60,0,255:48:68:63
#
#    then you would set --ref AL513382
#
# Or if your VCF file looks like this:
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	2332STDY5573284
#CT18_ref	744	.	A	G	0.213	.	DP=49;VDB=5.579080e-02;RPB=9.166667e-01;AF1=0.495;AC1=1;DP4=6,24,9,0;MQ=60;FQ=-10.4;PV4=2.4e-05,7.9e-08,1,1	GT:PL:DP:SP:GQ	0/1:17,0,255:39:46:5
#CT18_ref	762	.	A	G	29	.	DP=43;VDB=1.182296e-01;RPB=1.876799e+00;AF1=0.5;AC1=1;DP4=6,23,12,0;MQ=60;FQ=32;PV4=2.4e-06,1.4e-12,1,1	GT:PL:DP:SP:GQ	0/1:59,0,255:41:56:62
#CT18_ref	3727	.	A	G	30	.	DP=51;VDB=1.213276e-01;RPB=1.690455e+00;AF1=0.5;AC1=1;DP4=7,26,15,0;MQ=60;FQ=33;PV4=1.6e-07,2.4e-17,1,1	GT:PL:DP:SP:GQ	0/1:60,0,255:48:68:63
#
#    then you would set --ref CT18_ref
#
# Note CT18 has the genotype 3.2.1. If you specify the wrong reference name for the VCF
#   file, the script will fail to find any confident SNPs that differentiate your input 
#   strain from CT18 and will report that your strain is 3.2.1. So if all your strains
#   genotype as 3.2.1, check that you have 
#
# Authors - Kat Holt (kholt@unimelb.edu.au)
#
# Dependencies:
#	 NONE
#
# Last modified - Nov 9, 2015
#

from argparse import (ArgumentParser, FileType)
import os, sys, re, collections, operator
import gzip

def parse_args():
	"Parse the input arguments, use '-h' for help"
	parser = ArgumentParser(description='VCF to Typhi genotypes')
	parser.add_argument(
		'--vcf', nargs='+', type=str, required=False,
		help='VCF file(s) to genotype (Mapping MUST have been done using CT18 as a reference sequence)')
	parser.add_argument( '--ref', type=str, required=True, help='Name of the reference in the VCF file (#CHROM column). Note that CT18 has genotype 3.2.1. If all your strains return this genotype, it is likely you have specified the name of the refrence sequence incorrectly; please check your VCFs.') 
	parser.add_argument( '--phred', type=int, default=20, help='Minimum phred quality to count a variant call vs CT18 as a true SNP (default 20)')
	parser.add_argument( '--min_prop', type=float, default=0.1, help='Minimum proportion of reads required to call a SNP (default 0.1)')
	return parser.parse_args()


### SNP definitions

loci = [655112, 773487, 1804415, 1840727, 3640678, 270120, 102135, 316489, 4105384, 555826, 2360997, 4664137, 2166082, 30192, 4288272, 2737027, 1215983, 4132985, 146673, 2517324, 3920009, 3276735, 1173984, 2683312, 4013386, 4094437, 827036, 3476114, 4355243, 2723847, 4388609, 703762, 431216, 3095443, 316186, 1934711, 2811222, 3092900, 2723724, 3437570, 1780319, 1535365, 1792810, 3062270, 1799842, 432732, 3069182, 2732615, 3770391, 2269835, 4215341, 4602946, 3368641, 2245432, 3164162, 3923165, 1811809, 3729635, 3817752, 183033, 1615350, 2342045, 3996717, 2640029, 989024, 3806278, 1611156, 2348902]
snp_alleles = ['T', 'A', 'C', 'A', 'A', 'T', 'A', 'T', 'T', 'A', 'T', 'T', 'A', 'T', 'T', 'G', 'G', 'A', 'A', 'A', 'T', 'A', 'T', 'A', 'A', 'A', 'A', 'T', 'C', 'A', 'C', 'A', 'A', 'A', 'T', 'T', 'T', 'T', 'G', 'C', 'A', 'T', 'T', 'C', 'C', 'T', 'G', 'A', 'T', 'G', 'C', 'T', 'T', 'A', 'A', 'A', 'T', 'T', 'T', 'T', 'T', 'A', 'T', 'A', 'T', 'A', 'A', 'T']
groups = ['0.1', '0.0.1', '0.0.2', '0.0.3', '0.1.1', '0.1.2', '0.1.3', '1', '1.1', '1.1.1', '1.1.2', '1.1.3', '1.1.4', '1.2', '1.2.1', '2', '2.0.1', '2.0.2', '2.1', '2.1.1', '2.1.2', '2.1.3', '2.1.4', '2.1.5', '2.1.6', '2.1.7', '2.1.8', '2.1.9', '2.2', '2.2.1', '2.2.2', '2.2.3', '2.2.4', '2.3', '2.3.1', '2.3.2', '2.3.3', '2.3.4', '2.3.5', '2.4', '2.4.1', '2.5', '2.5.1', '3', '3.0.1', '3.0.2', '3.1', '3.1.1', '3.1.2', '3.2', '3.2.1', '3.2.2', '3.3', '3.3.1', '3.4', '3.5', '3.5.1', '3.5.2', '3.5.3', '3.5.4', '4', '4.1', '4.1.1', '4.2', '4.2.1', '4.2.2', '4.2.3', '4.3']

# check if this SNP defines a group
def checkSNP(vcf_line_split, this_groups, proportions, args):
	snp = int(vcf_line_split[1])
	if snp in loci:
		i = loci.index(snp)
		if float(vcf_line_split[5]) > args.phred:
			m = re.search("DP4=(\d+),(\d+),(\d+),(\d+)",vcf_line_split[7])
			if m != None:
				alt_read_count = int(m.group(3)) + int(m.group(4))
				total_read_count = alt_read_count + int(m.group(1)) + int(m.group(2))
				snp_proportion = float(alt_read_count)/total_read_count
			else:
				snp_proportion = float(-1) # set unknowns to negative so that we know this is not a real proportion
			if snp_proportion > args.min_prop:
				this_allele = vcf_line_split[4]
				if this_allele == snp_alleles[i]:
					this_groups.append(groups[i]) # retrieve the group that this SNP defines
					proportions[groups[i]] = snp_proportion
	return (this_groups, proportions)

# sort groups into the three levels (primary, clade, subclade)
def parseGeno(this_groups, proportions):
	subclades = []
	clades = []
	primary = []
	for group in this_groups:
		level = len(group.split("."))
		if level == 3:
			subclades.append(group)
		elif level == 2:
			clades.append(group)
		elif level == 1:
			primary.append(group)
			
	# fix 2.3, 2.2 nesting
	if ('2.2' in clades) and ('2.3' in clades):
		clades.remove('2.2')
				
	# fix 3.5.3, 3.5.4 nesting
	if ('3.5.3' in subclades) and ('3.5.4' in subclades):
		subclades.remove('3.5.3')
		
	# fix 2.3.1, 2.3.3 nesting
	if ('2.3.1' in subclades) and ('2.3.2' in subclades):
		subclades.remove('2.3.2')

	# fix 2.3.1, 2.3.3 nesting
	if ('2.3.5' in subclades) and ('2.3.3' in subclades):
		subclades.remove('2.3.3')
				
	# fix primary clades relative to CT18 = 3.2.1, ie has clade1, clade2, clade3 SNPs
	if len(primary) == 1:
		if '3' in primary:
			primary = ['2'] # clade 2 differs from CT18 by the clade3-defining SNP
			# note other option is clade 4 snp, which defines clade 4 relative to CT18
	elif len(primary) == 2:
		if ('2' in primary) and ('3' in primary):
			primary = ['1'] # clade 2 differs from CT18 by the clade3-defining SNP
	elif len(primary) == 0:
		primary = ['3']
	elif len(primary) == 3:
		if ('1' in primary) and ('2' in primary) and ('3' in primary):
			primary = ['0']
	
	# fix clade relative to CT18:
	if '3.2' in clades:
		clades.remove('3.2') # anything NOT in 3.2 will have this SNP
	else:
		if len(clades) == 0:
			clades.append('3.2') # anything with no clade, and 3.2 SNP not called, belongs in 3.2 with CT18

	# fix 3.5.3, 3.5.4 nesting
	if ('3.5.3' in clades) and ('3.5.4' in clades):
		clades.remove('3.5.3')
				
	# fix subclades relative to CT18:
	if '3.2.1' in subclades:
		subclades.remove('3.2.1') # anything NOT in 3.2.1 will have this SNP
	else:
		if len(subclades) == 0:
			subclades.append('3.2.1') # anything with no subclade, and 3.2.1 SNP NOT called, belongs in 3.2.1 with CT18
		
	# add zero-th clade/subclade where unresolved -- disabled
	#if len(clades) == 0:
	#	if len(primary) == 1:
	#		clades.append(primary[0] + ".0")
	#if len(subclades) == 0:
	#	if len(clades) == 1:
	#		subclades.append(clades[0] + ".0")
			
	info = ",".join(subclades) + "\t" + ",".join(clades) + "\t" + ",".join(primary)
	
	# add proportion of reads supporting each of these groups
	p = []
	for group in subclades:
		if group in proportions:
			p.append(str(round(proportions[group],2)))
	info += "\t" + ",".join(p)
	
	p = []
	for group in clades:
		if group in proportions:
			p.append(str(round(proportions[group],2)))
	info += "\t" + ",".join(p)
	
	p = []
	for group in primary:
		if group in proportions:
			p.append(str(round(proportions[group],2)))
	info += "\t" + ",".join(p)
	
	return info

# main function
def main():
	args = parse_args()

	print "\t".join(["File","Subclade","Clade","PrimaryClade","Support_Subclade","Support_Clade","Support_PrimaryClade"])
	
	for vcf in args.vcf:
		this_groups = [] # list of groups identified by defining SNPs
		proportions = {} # proportion of reads supporting each defining SNP; key = group, value = proportion
		
		# read file
		(file_name,ext) = os.path.splitext(vcf)
		if ext == "gz":
			f = gzip.open(file_name,"r")
		else:
			f = open(vcf,"r")

		any_ref_line = 0
		
		for line in f:
			if not line.startswith("#"):
				x = line.rstrip().split()
				if x[0] == args.ref:
					# parse this SNP line
					(this_groups, proportions) = checkSNP(x, this_groups, proportions, args)
		
		f.close()

		if any_ref_line > 0:
			info = parseGeno(this_groups, proportions)
			print vcf + "\t" + info
		else:
			print vcf + "\twrongref"
				

# call main function
if __name__ == '__main__':
	main()
