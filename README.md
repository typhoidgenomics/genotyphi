# genotyphi

This script assigns genotypes to Salmonella Typhi genomes. The genotyping scheme is described in this paper (please cite if you use this method): ["An extended genotyping framework for Salmonella enterica serovar Typhi, the cause of human typhoid", Wong et al, 2016, Nature Communications](http://www.nature.com/articles/ncomms12827/). It is also summarised in [this blog post](https://holtlab.net/2016/10/12/global-picture-typhoid/).

Inputs are BAM or VCF files (mapped to the Typhi CT18 reference genome, [AL513382.1](https://www.ncbi.nlm.nih.gov/nuccore/AL513382.1)).

For short read data, we recommend using the raw alignments (BAM files) as input (via --bam), since VCFs generated by different pipelines may have filtered out important data at the genotyping loci. Note you also need to supply the reference sequence that was used. If you don't have BAMs or you are really confident in the quality of your SNP calls in existing VCF files, you can provide your own VCFs (generated by read mapping) via --vcf.

For assemblies, we recommend using [ParSNP](http://harvest.readthedocs.org/) (version 1.0.1) to align genomes to the CT18 reference. The resulting multi-sample VCF file(s) can be passed directly to this script via --vcf_parsnp.

Dependencies: Python 2.7.5+ ([SAMtools](http://samtools.sourceforge.net/) and [BCFtools](https://samtools.github.io/bcftools/) are also required if you are working from BAM files.  Genotyphi has been tested with versions 1.1, 1.2, and 1.3 of both SAMtools and bcftools, subsequently we advise using the same version of both of these dependencies together i.e. SAMtools v1.2 and bcftools v1.2).

[Basic Usage](https://github.com/katholt/genotyphi/#basic-usage---own-bam-recommended-if-you-have-reads)

[Options](https://github.com/katholt/genotyphi/#options)

[Outputs](https://github.com/katholt/genotyphi/#outputs)

[Generating input BAMS from reads (with example)](https://github.com/katholt/genotyphi/#generating-input-bams-from-reads)

[Generating input VCFs from assemblies (with example)](https://github.com/katholt/genotyphi/#generating-input-vcfs-from-assemblies)

## Basic Usage - own BAM (recommended if you have reads)

Note the BAM files must be sorted (e.g. using samtools sort)

```
python genotyphi.py --mode bam --bam *.bam --ref AL513382.fasta --ref_id AL513382.1 --output genotypes.txt
```

## Basic Usage - own VCF

```
python genotyphi.py --mode vcf --vcf *.vcf --ref_id AL513382 --output genotypes.txt
```

## Basic Usage - assemblies aligned with ParSNP (recommended if you only have assembly data available and no reads)

```
python genotyphi.py --mode vcf_parsnp --vcf parsnp.vcf --output genotypes.txt
```

## Options

### Required options

```
-- mode            Run mode, either bam, vcf or vcf_parsnp
```

### Mode specific options

#### --mode bam

Requires [SAMtools](http://samtools.sourceforge.net/) and [BCFtools](https://samtools.github.io/bcftools/)

```
--bam              Path to one or more BAM files, generated by mapping reads to the Typhi CT18 reference genome (AL513382)
                   Note the SNP coordinates used here for genotyping are relative to Typhi CT18 (AL513382) 
                   so the input MUST be a BAM obtained via mapping to this reference sequence.

--ref              Reference sequence file used for mapping (fasta).

--ref_id           Name of the Typhi CT18 chromosome reference in your VCF file.
                   This is the entry in the first column (#CHROM) of the data part of the file.
                   This is necessary in case you have mapped to multiple replicons (e.g. chromosome and
                   plasmid) and all the results appear in the same VCF file.

--samtools_location     Specify the location of the folder containing the SAMtools installation if not standard/in path [optional]

--bcftools_location     Specify the location of the folder containing the bcftools installation if not standard/in path [optional]
```

#### --mode vcf

```
--vcf              Path to one or more VCF files, generated by mapping reads to the Typhi CT18 reference genome (AL513382)
                   Note the SNP coordinates used here for genotyping are relative to Typhi CT18 (AL513382) 
                   so the input MUST be a VCF obtained via mapping to this reference sequence.
```

#### --mode vcf_parsnp

```
--vcf              Path to one or more VCF files generated by mapping assemblies to the Typhi CT18 reference genome (AL513382)
                   using ParSNP (--ref_id is optional, default value is '1').
```

### Other options

```
--phred                 Minimum phred quality to count a variant call vs CT18 as a true SNP (default 20)

--min_prop              Minimum proportion of reads required to call a SNP (default 0.1)

--output                Specify the location of the output text file (default genotypes_[timestamp].txt)

```

## Outputs

Output is to standard out, in tab delimited format.

The script works by looking for SNPs that define three nested levels of phylogenetic lineages for S. Typhi: primary clusters, clades and subclades. These are named in the form 1.2.3, which indicates primary cluster 1, clade 1.2, and subclade 1.2.3.

All genotypes implicated by the detected SNPs will be reported (comma separated if multiple are found). 

Usually subclade, clade and primary cluster SNPs will be compatible, i.e. you would expect to see 1.2.3 in the subclade column, 1.2 in the clade column and 1 in the primary clade column.

The first column 'Final_call' indicates the highest resolution genotype (subclade > clade > primary clade) that could be called.

As recombination is extremely rare in S. Typhi, it is unlikely that DNA isolated from a single S. Typhi culture would have high quality SNPs that are compatible with multiple genotypes, or incompatible clade/subclade combinations... therefore if you see this in the output, it is worth investigating further whether you have contaminated sequence data or a genuine recombination. To provide very basic assistance with this the output file indicates, for each genotype-defining SNP detected, the proportion of reads that support the genotype call. A genuine recombinant strain would have p=1 for multiple incompatible SNPs, whereas a mixed DNA sample might have p~0.5. Note that if you are genotyping from assemblies, we have no information on read support to work with; instead you will see an 'A' in the Final_call_support column to indicate this result comes from a genome assembly.

WARNING: Note the reference genome CT18 has the genotype 3.2.1. It is therefore possible that unexpected behaviour (e.g. wrong reference, unspecified reference name, or problems that results in overall poor quality SNP calls in the VCF) can sometimes result in calls of 3.2.1 ... so if you see lots of strains being assigned this genotype, it is worth investigating further.

## Generating input BAMS from reads

We recommend using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to align reads to the CT18 reference, and [SAMtools](http://http://samtools.sourceforge.net/) to convert the *.sam file to the *.bam format.  The resulting bam file(s) can be passed directly to this script via --bam.  Please note the differences in the commands listed below when using SAMtools v1.1/1.2 vs. SAMtools v1.3.1 to generate bam files.

For example:

```
# Download CT18 (AL513382) reference genome

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000195995.1_ASM19599v1/GCA_000195995.1_ASM19599v1_genomic.fna.gz

gunzip GCA_000195995.1_ASM19599v1_genomic.fna.gz
mv GCA_000195995.1_ASM19599v1_genomic.fna CT18.fasta

# Download reads for S. Typhi 8(04)N, isolated from Vietnam in 2004

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR343/ERR343343/ERR343343_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR343/ERR343343/ERR343343_2.fastq.gz

# Use bowtie to map reads to the CT18 reference genome

For example, to align paired end reads to the CT18 reference genome sequence:

bowtie2-build CT18.fasta CT18

bowtie2 -p 2 -x CT18 -1 ERR343343_1.fastq.gz -2 ERR343343_2.fastq.gz -S output.sam
 
samtools view -bS output.sam > unsorted_output.bam

samtools sort unsorted_output.bam output

(or, 'samtools sort unsorted_output.bam > output.bam' for SAMtools v1.3.1 instead of SAMtools v1.2/1.1)

# Call Typhi genotypes from the resulting BAM(s)

python genotyphi.py --mode bam --bam output.bam --ref CT18.fasta --ref_id AL513382.1 --output genotypes_test.txt

```

#### Output

For bam output the column 'Number of SNPs called' reports the number of SNPs present in the resultant VCF file(s).

```
File    Final_call  Final_call_support  Subclade    Clade   PrimaryClade    Support_Subclade    Support_Clade   Support_PrimaryClade    Number of SNPs called
output.vcf      4.3     0.97            4.3     4               0.97    1.0     68

```

## Generating input VCFs from assemblies

We recommend using [ParSNP](http://harvest.readthedocs.org/) to align genomes to the CT18 reference. The resulting multi-sample VCF file(s) can be passed directly to this script via --vcf_parsnp.

For example:

```
# Download CT18 reference genome

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000195995.1_ASM19599v1/GCA_000195995.1_ASM19599v1_genomic.fna.gz

gunzip GCA_000195995.1_ASM19599v1_genomic.fna.gz
mv GCA_000195995.1_ASM19599v1_genomic.fna CT18.fasta

# Download two example Typhi genomes for genotyping

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Salmonella_enterica/latest_assembly_versions/GCA_000245535.1_ASM24553v1/GCA_000245535.1_ASM24553v1_genomic.fna.gz
wget ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/Salmonella_enterica/latest_assembly_versions/GCA_000007545.1_ASM754v1/GCA_000007545.1_ASM754v1_genomic.fna.gz

gunzip GCA_000245535.1_ASM24553v1_genomic.fna.gz
gunzip GCA_000007545.1_ASM754v1_genomic.fna.gz

mkdir genomes/
mv GCA_000245535.1_ASM24553v1_genomic.fna genomes/Pstx12.fasta
mv GCA_000007545.1_ASM754v1_genomic.fna genomes/Ty2.fasta

# Use ParSNP to generate variant calls (VCF) for these genomes against the CT18 reference sequence

parsnp -r CT18.fasta -d genomes/ -o output

# Call Typhi genotypes from the resulting VCF

python genotyphi.py --mode vcf_parsnp --vcf output/parsnp.vcf --output genotypes_parsnptest.txt

```
#### Output

P-stx-12 should be called as clade 4.3 (otherwise known as H58), Ty2 should be called as clade 4.1. 

The 'A' in the Final_call_support indicates the genotype calls were made from assembled data, hence no read support information is available.

```
File	Final_call	Final_call_support	Subclade	Clade	PrimaryClade	Support_Subclade	Support_Clade	Support_PrimaryClade
Pstx12.fasta	4.3	A		4.3	4			
Ty2.fasta	4.1	A		4.1	4	

```
