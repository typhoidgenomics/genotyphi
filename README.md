# genotyphi

The GenoTyphi genotyping scheme divides the *Salmonella* Typhi population into 4 major lineages, and >75 different clades and subclades. The scheme specification is detailed in the file 'Genotype_specification.csv' in this repository, which also includes the standard clade-level colour codes that we use for consistency across papers.

The scheme is described in this paper, **which serves as the primary citation for the scheme**: ["An extended genotyping framework for Salmonella enterica serovar Typhi, the cause of human typhoid", Wong et al, 2016, Nature Communications](http://www.nature.com/articles/ncomms12827/).

Subsequent updates to the genotyping scheme, including new genotypes and mutations conferring resistance to fluoroquinolones and azithromycin, are summarised in ["Five years of GenoTyphi: updates to the global Salmonella Typhi genotyping framework", Dyson & Holt, 2021, Journal of Infectious Diseases](https://doi.org/10.1093/infdis/jiab414).

Two implementations of the code that can assign genotypes to Typhi genomes are available in this repo:

* [Mykrobe implementation](https://github.com/katholt/genotyphi/blob/main/README.md#mykrobe-implementation), which takes as input **sequence reads (fastq files)** from Illumina or long-read platforms. It also detects known mutations in the quinolone-resistance determining region (QRDR) of genes *gyrA*, *gyrB* and *parC*; and the *acrB*-R717Q/L mutations associated with azithromycin resistance (AMR SNVs are listed in 'AMR_mutation_alleles.csv' in this repository); acquired antimicrobial resistance genes; plasmid replicons and major subtypes of the IncHI1 plasmid typically associated with multidrug resistance (list of targets is in 'AMR_genes_drugs.txt').

* [Original Python implementation](https://github.com/katholt/genotyphi/blob/main/README.md#original-implementation-pre-mapped-data), which takes as input **BAM or VCF files that the user has already generated** by mapping Illumina reads to the reference genome CT18. It also detects the QRDR and *acrB* mutations listed in 'AMR_mutation_alleles.csv' in this repository.

The GenoTyphi scheme is also available via the online analysis platform [Pathogenwatch](https://pathogen.watch/), which facilitates automated analysis of Typhi genomic reads or assemblies as described in [this paper](https://www.nature.com/articles/s41467-021-23091-2).

**Whichever tool you use to access the GenoTyphi scheme, please cite the [GenoTyphi paper](http://www.nature.com/articles/ncomms12827/).**

**If you use the scripts in this repository, please also cite the repository:** [![DOI](https://zenodo.org/badge/45819844.svg)](https://zenodo.org/badge/latestdoi/45819844)

## Mykrobe Implementation

### Installing mykrobe

First, install Mykrobe (v0.10.0+) as per the instructions on the [Mykrobe github](https://github.com/Mykrobe-tools/mykrobe).

Once Mykrobe is installed, you can run the following two commands to ensure you have the most up-to-date panels for genotyping, including the [Typhi panel](https://doi.org/10.26180/14881701.v1):
```
mykrobe panels update_metadata
mykrobe panels update_species all
```

You can check what version of the scheme is currently loaded in your Mykrobe installation via:
```
mykrobe panels describe
```

While the above description will guarantee up-to-date versions of both mykrobe and the genotyphi panel for mykrobe, it may be difficult or confusing if you are not familiar with docker. 

If you are not familiar with docker, you can install a version of mykrobe which is pre-configured for Typhi by running:

```
docker pull flashton/mykrobe_for_typhi:latest
```

The version of mykrobe inside this container is `mykrobe:0.10.0--py39h98c8e45_0`, and the genotyphi catalogue it uses was downloaded on 2021-11-11.

### Running mykrobe

Inputs are fastq files.

Mykrobe can be run on each individual sample using the command below. Replace `aSample` with the name of your isolate. The command below uses Illumina data.

```
mykrobe predict --sample aSample --species typhi --format json --out aSample.json --seq aSample_1.fastq.gz aSample_2.fastq.gz
```

To run on ONT data instead, add the `--ont` flag to your command.

Further details on options can be found on the Mykrobe wiki: https://github.com/Mykrobe-tools/mykrobe/wiki

If you are running the version of mykrobe from `flashton/mykrobe_for_typhi:latest`, you can run mykrobe on your Typhi data with a single command. For example:

```
docker run --rm -v /path/to/your/data:/data flashton/mykrobe_for_typhi:latest mykrobe predict --sample aSample --species typhi --format json --out /data/aSample.json --seq /data/my_fastq.gz --ont
```

Note that `my_fastq.gz` has to be inside the `/path/to/your/data` directory for this command to work.


### Parse mykrobe output

#### Code

`parse_typhi_mykrobe.py`

We have provided a custom python3 script, parse_typhi_mykrobe.py, that will take a group of JSON files output by Mykrobe and summarise these into a single, tab-delimited table.

The parser script will only report details of calls for genomes that are identified by Mykrobe as Typhi. Currently, a sample must have >=85% identity to Typhi MLST locus sequences to be called by the parser. This threshold may not be low enough to correctly parse json files created by analysing ONT data (however all the Mykrobe calls will still be present in the json file).

Note that due to the nested hierarchical nature of the GenoTyphi scheme, we needed to create fake levels within the scheme to facilitate correct calling by Mykrobe. These fake levels are not reported in the output generated by the parser, but they are present in the raw json files output by Mykrobe. These can be identified in the json output as they are prepended by the word "lineage", and will always have a call of 0 from Mykrobe.

#### Dependencies

* Python 3 
* library `pandas`

#### Input

* json files output by Mykrobe
* prefix for output file

#### Output

* tsv file, one row per json file

#### Example command
```
python parse_typhi_mykrobe.py --jsons *.json --prefix mykrobe_out
```

### Explanation of columns in the output:
* **genome**: sample ID
* **species**: species call from Mykrobe (Typhi or unknown; determined by matching to Typhi STs from the 7-locus MLST scheme)
* **spp_percent:** percentage coverage to the Typhi ST probes
* **final genotype**: final genotype call from Mykrobe, using the latest version of GenoTyphi
* **confidence**: measure of confidence in the final genotype call, summarising read support across all levels in the hierarchy (lineage, clade, subclade, etc)
  * _strong_ - high quality calls (quality of '1' reported by Mykrobe) for ALL levels;
  * _moderate_ - reduced confidence for ONE node (Mykrobe quality of '0.5', but still with >50% read support for the marker allele), high quality calls for ALL OTHERS;
  * _weak_ - low quality for one or more nodes (Mykrobe quality of '0.5' and <50% read support OR Mykrobe quality of '0').
* **acrB_R717L, acrB_7171Q**: calls for individual *acrB* (azithromycin resistance) mutation. 0 indicates mutation is absent, 1 indicates mutation is present.
* **num QRDR**: Total number of mutations detected in the quinolone-resistance determining regions (QRDR) of genes _gyrA_, _parC_ and _gyrB_.
* **lowest support for genotype marker**: For any markers in the final genotype call that do not have a Mykrobe quality of '1', this column reports the percentage of reads supporting the marker allele at the most poorly supported marker (details of all such markers appear in the 'poorly supported markers column').
* **poorly supported markers**: Lists any markers in the final genotype call that do not have Mykrobe quality of '1'. Markers are separated by ';', values in brackets represent the quality call from Mykrobe, followed by the read depth at the alternate / reference alleles. The lowest read support amongst these markers is reported in the previous column.
* **max support for additional markers**: For any markers detected that are incongruent with the final genotype call, this column reports the percentage of reads supporting the marker allele at the best supported additional marker.
* **additional markers**: Lists any markers that are incongruent with the final genotype call. Markers are separated by ';', and the format is identical to column _poorly supported markers_. The highest read support for any such marker is reported in the previous column.
* **node support**: A list of all markers in the final genotype call with their Mykrobe quality calls (1, 0.5, or 0) and the read depths at the marker allele / reference allele (as per _poorly supported markers_ and _additional markers_).
* **parC_S80R, parC_S80I, parC_E84G, parC_E84K, gyrA_S83F, gyrA_S83Y, gyrA_D87G, gyrA_D87N, gyrA_D87V, gyrA_D87Y, gyrB_S464F, gyrB_S464Y**: calls for each individual QRDR mutation. 0 indicates mutation is absent, 1 indicates mutation is present.
* **remaining columns** indicate presence (1) or absence (0) of each AMR gene or plasmid replicon as indicated by the header of the column. For the AMR genes, _specific alleles are **not detected**_ (with the exception of dfrA, where Mykrobe can distinguish between dfrA7, dfrA5 and dfrA15).
  * **IncHI1_ST6** indicates whether the IncHI1 replicon (from column IncHI1BR27) is pST6 (1) or not (0).

Sequences and details of probes are available [here](https://doi.org/10.26180/14881701.v1).


## Original Implementation (pre-mapped data)

### Code

`genotyphi.py`

### Dependencies
Python 2.7.5+

[SAMtools](http://samtools.sourceforge.net/) and [BCFtools](https://samtools.github.io/bcftools/) are also required if you are working from BAM files.  Genotyphi has been tested with versions 1.1, 1.2, and 1.3 of both SAMtools and bcftools, subsequently we advise using the same version of both of these dependencies together i.e. SAMtools v1.2 and bcftools v1.2

### Input files
Inputs are BAM or VCF files (mapped to the Typhi CT18 reference genome, [AL513382.1](https://www.ncbi.nlm.nih.gov/nuccore/AL513382.1)).

For short read data, we recommend using the raw alignments (BAM files) as input (via --bam), since VCFs generated by different pipelines may have filtered out important data at the genotyping loci. Note you also need to supply the reference sequence that was used. If you don't have BAMs or you are really confident in the quality of your SNP calls in existing VCF files, you can provide your own VCFs (generated by read mapping) via --vcf.

For assemblies, we recommend using [ParSNP] and [HarvestTools] (http://harvest.readthedocs.org/) (version 1.0.1) to align genomes to the CT18 reference. The resulting multi-sample VCF file(s) can be passed directly to this script via --vcf_parsnp.

### Usage

[Basic Usage](https://github.com/katholt/genotyphi/#basic-usage---own-bam-recommended-if-you-have-reads)

[Options](https://github.com/katholt/genotyphi/#options)

[Outputs](https://github.com/katholt/genotyphi/#outputs)

[Generating input BAMS from reads (with example)](https://github.com/katholt/genotyphi/#generating-input-bams-from-reads)

[Generating input VCFs from assemblies (with example)](https://github.com/katholt/genotyphi/#generating-input-vcfs-from-assemblies)

#### Basic Usage - own BAM (recommended if you have reads)

Note the BAM files must be sorted (e.g. using samtools sort)

```
python genotyphi.py --mode bam --bam *.bam --ref AL513382.fasta --ref_id AL513382.1 --output genotypes.txt
```

#### Basic Usage - own VCF

```
python genotyphi.py --mode vcf --vcf *.vcf --ref_id AL513382 --output genotypes.txt
```

#### Basic Usage - assemblies aligned with ParSNP (recommended if you only have assembly data available and no reads)

```
python genotyphi.py --mode vcf_parsnp --vcf parsnp.vcf --output genotypes.txt
```

#### Options

##### Required options

```
-- mode            Run mode, either bam, vcf or vcf_parsnp
```

##### Mode specific options

##### --mode bam

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

##### --mode vcf

```
--vcf              Path to one or more VCF files, generated by mapping reads to the Typhi CT18 reference genome (AL513382)
                   Note the SNP coordinates used here for genotyping are relative to Typhi CT18 (AL513382) 
                   so the input MUST be a VCF obtained via mapping to this reference sequence.
```

##### --mode vcf_parsnp

```
--vcf              Path to one or more VCF files generated by mapping assemblies to the Typhi CT18 reference genome (AL513382)
                   using ParSNP (--ref_id is optional, default value is '1').
```

#### Other options

```
--phred                 Minimum phred quality to count a variant call vs CT18 as a true SNP (default 20)

--min_prop              Minimum proportion of reads required to call a SNP (default 0.1)

--output                Specify the location of the output text file (default genotypes_[timestamp].txt)

```

#### Outputs

Output is to standard out, in tab delimited format.

The script works by looking for SNPs that define nested levels of phylogenetic lineages for S. Typhi: primary clusters, clades and subclades. These are named in the form 1.2.3, which indicates primary cluster 1, clade 1.2, and subclade 1.2.3.  For H58 the script will call lineages 1 and 2 as 4.3.1.1 and 4.3.1.2.

All genotypes implicated by the detected SNPs will be reported (comma separated if multiple are found). 

Usually subclade, clade and primary cluster SNPs will be compatible, i.e. you would expect to see 1.2.3 in the subclade column, 1.2 in the clade column and 1 in the primary clade column.

The first column 'Final_call' indicates the highest resolution genotype (subclade > clade > primary clade) that could be called.

As recombination is extremely rare in S. Typhi, it is unlikely that DNA isolated from a single S. Typhi culture would have high quality SNPs that are compatible with multiple genotypes, or incompatible clade/subclade combinations... therefore if you see this in the output, it is worth investigating further whether you have contaminated sequence data or a genuine recombination. To provide very basic assistance with this the output file indicates, for each genotype-defining SNP detected, the proportion of reads that support the genotype call. A genuine recombinant strain would have p=1 for multiple incompatible SNPs, whereas a mixed DNA sample might have p~0.5. Note that if you are genotyping from assemblies, we have no information on read support to work with; instead you will see an 'A' in the Final_call_support column to indicate this result comes from a genome assembly.

WARNING: Note the reference genome CT18 has the genotype 3.2.1. It is therefore possible that unexpected behaviour (e.g. wrong reference, unspecified reference name, or problems that results in overall poor quality SNP calls in the VCF) can sometimes result in calls of 3.2.1 ... so if you see lots of strains being assigned this genotype, it is worth investigating further.

#### Generating input BAMS from reads

We recommend using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to align reads to the CT18 reference, and [SAMtools](http://http://samtools.sourceforge.net/) to convert the *.sam file to the *.bam format.  The resulting bam file(s) can be passed directly to this script via --bam.  Please note the differences in the commands listed below when using SAMtools v1.1/1.2 vs. SAMtools v1.3.1 to generate bam files.

For example:

```
# Download CT18 (AL513382) and unzip the CT18 reference genome

wget -O CT18.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Salmonella_enterica/reference/GCA_000195995.1_ASM19599v1/GCA_000195995.1_ASM19599v1_genomic.fna.gz

gunzip CT18.fasta.gz

# Separate the chromosome sequence from the plasmids with the emboss toolkit
seqretsplit CT18.fasta

mv nc_003198.1.fasta CT18.fasta

# Replace the header line of the CT18.fasta file with the reference id i.e. chage

>NC_003198.1 Salmonella enterica subsp. enterica serovar Typhi str. CT18, complete genome

# to
> AL513382.1

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

Expected output:

For bam output the column 'Number of SNPs called' reports the number of SNPs present in the resultant VCF file(s).

```
File	Final_call	Final_call_support	Subclade	Clade	PrimaryClade	Support_Subclade	Support_Clade	Support_PrimaryClade	Number of SNPs called	QRDR mutations
output.vcf	4.3.1.2	1.0	4.3.1.2		4	1.0		1.0	77	 gyrA-D87G

```

#### Generating input VCFs from assemblies

We recommend using [ParSNP](http://harvest.readthedocs.org/) to align genomes to the CT18 reference. The resulting multi-sample VCF file(s) can be passed directly to this script via --vcf_parsnp.

For example:

```
# Download CT18 reference genome

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Salmonella_enterica/reference/GCA_000195995.1_ASM19599v1/GCA_000195995.1_ASM19599v1_genomic.fna.gz

gunzip GCA_000195995.1_ASM19599v1_genomic.fna.gz

# Separate CT18 chromosome with emboss
seqretsplit GCA_000195995.1_ASM19599v1_genomic.fna
mv al513382.1.fasta CT18.fasta

# Download two example Typhi genomes for genotyping

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/245/535/GCA_000245535.1_ASM24553v1/GCA_000245535.1_ASM24553v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/007/545/GCA_000007545.1_ASM754v1/GCA_000007545.1_ASM754v1_genomic.fna.gz

gunzip GCA_000245535.1_ASM24553v1_genomic.fna.gz
gunzip GCA_000007545.1_ASM754v1_genomic.fna.gz

mkdir genomes/
mv GCA_000245535.1_ASM24553v1_genomic.fna genomes/Pstx12.fasta
mv GCA_000007545.1_ASM754v1_genomic.fna genomes/Ty2.fasta

# Use ParSNP and Harvest tools to generate variant calls (VCF) for these genomes against the CT18 reference sequence

parsnp -r CT18.fasta -d genomes/ -o output
harvesttools -i ./output/parsnp.ggr -V parsnp.vcf

# Call Typhi genotypes from the resulting VCF

python genotyphi.py --mode vcf_parsnp --vcf parsnp.vcf --ref_id AL513382.1 --output genotypes_parsnptest.txt

```
Expected Output:

P-stx-12 should be called as clade 4.3.1 (otherwise known as H58), Ty2 should be called as clade 4.1. 

The 'A' in the Final_call_support indicates the genotype calls were made from assembled data, hence no read support information is available.

```
File	Final_call	Final_call_support	Subclade	Clade	PrimaryClade	Support_Subclade	Support_Clade	Support_PrimaryClade
Pstx12.fasta	4.3.1	A		4.3.1	4			
Ty2.fasta	4.1	A		4.1	4	

```
