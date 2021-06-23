# genotyphi

#### Intro & citation
This script assigns genotypes to *Salmonella* Typhi genomes and detectes known mutations in the quinolone-resistance determining region (QRDR) of genes *gyrA* and *parC*, as well as the acrB-R717Q/L mutations confering resistance to azithromycin. 

The genotyping scheme is described in this paper, **which serves as the primary citation for the scheme**: ["An extended genotyping framework for Salmonella enterica serovar Typhi, the cause of human typhoid", Wong et al, 2016, Nature Communications](http://www.nature.com/articles/ncomms12827/).

The orignal genotyping method is also summarised in [this blog post](https://holtlab.net/2016/10/12/global-picture-typhoid/).

Subsequent updates to the genotyping scheme are summarised in [this preprint](https://www.biorxiv.org/content/10.1101/2021.04.28.441766v1) and below.

[![DOI](https://zenodo.org/badge/45819844.svg)](https://zenodo.org/badge/latestdoi/45819844)

#### Additional features & citations

An updated version of the code that **detects QRDR mutations, and genotypes H58 lineages 4.3.1.1, 4.3.1.2 and 4.3.1.3** and sublineages 4.3.1.1.EA1, 4.3.1.2.EA2, 4.3.1.2.EA3, and 4.3.1.3.Bdq was developed and described in the papers ["Laboratory and Molecular Surveillance of Paediatric Typhoidal Salmonella in Nepal: Antimicrobial Resistance and Implications for Vaccine Policy"](https://journals.plos.org/plosntds/article?rev=1&id=10.1371/journal.pntd.0006408), Britto et al 2018, PLoS NTDs ["Population structure and antimicrobial resistance patterns of Salmonella Typhi isolates in Bangladesh from 2004 to 2016"](https://journals.plos.org/plosntds/article?rev=1&id=10.1371/journal.pntd.0006408), Rahman et al 2020, PLoS NTDs, and ["Multiple introductions of multidrug-resistant typhoid associated with acute infection and asymptomatic carriage, Kenya"](https://www.biorxiv.org/content/10.1101/2021.03.10.434750v1), biorxiv. Please cite all these in addition to the [primary paper](http://www.nature.com/articles/ncomms12827/) if you report genotyphi-inferred QRDR mutations and/or H58 sublineages 1,2,3 (ie 4.3.1.1, 4.3.1.2, 4.3.1.3) and East African H58 sublineages (4.3.1.1.EA1, 4.3.1.2.EA2, and 4.3.1.2.EA3). 

Following the publication of ["Emergence of an Extensively Drug-Resistant Salmonella enterica Serovar Typhi Clone Harboring a Promiscuous Plasmid Encoding Resistance to Fluoroquinolones and Third-Generation Cephalosporins", Klemm et al 2018, MBio](https://mbio.asm.org/content/9/1/e00105-18) genotyphi has been updated to identify members of the Pakistan XDR Typhi outbreak as 4.3.1.1.P1. Please cite this paper if you identify and report the presence of **genotype 4.3.1.1.P1**.

The GenoTyphi scheme is also available via the online analysis platform [PathogenWatch](https://pathogen.watch/), which facilitates automated analysis of Typhi genomic reads or assemblies as described in [this paper](https://www.biorxiv.org/content/10.1101/2020.07.03.186692v1.full). Whichever tool you use to access the GenoTyphi scheme, please do remember to cite the [GenoTyphi paper](http://www.nature.com/articles/ncomms12827/).

## Mykrobe Implementation: Inputs are raw reads

### Installing mykrobe

To use the Typhi panel, you need to install mykrobe from source (using pip), as the Typhi panel hasn't been incorporated into a full release yet (and so the conda instructions on the mykrobe installation page won't work). In the example below I have installed from source within a conda environment - if you don't want to create a conda environment, just skip to the third line.

```
conda create mykrobe
conda activate mykrobe
git clone https://github.com/Mykrobe-tools/mykrobe
cd mykrobe
pip3 install . && mykrobe panels update_metadata && mykrobe panels update_species all
```

You may run into issues with the dependency `mccortex`, please see the documentation on the Mykrobe wiki for more detailed instructions https://github.com/Mykrobe-tools/mykrobe/wiki/Installation#from-source

### Running mykrobe

Mykrobe can be run on each individual sample using the command below. Replace `aSample` with the name of your isolate. The command below uses Illumina data.

```
mykrobe predict --sample aSample --species typhi --format json --out aSample.json --seq aSample_1.fastq.gz aSample_2.fastq.gz
```

To run on ONT data instead, add the `--ont` flag to your command.

Further details on options can be found on the Mykrobe wiki: https://github.com/Mykrobe-tools/mykrobe/wiki

### Parse mykrobe output - script parse_typhi_mykrobe.py

We have provided a custom python3 script, parse_typhi_mykrobe.py, that will take a group of JSON files output by Mykrobe and summarise these into a single, tab-delimited table.

The parser script will only report details of calls for genomes that are identified by Mykrobe as Typhi. Currently, a sample must have >=85% identity to Typhi MLST locus sequences to be called by the parser. This threshold may not be low enough to correctly parse json files created by analysing ONT data (however all the Mykrobe calls will still be present in the json file).

Note that due to the nested hierarchical nature of the GenoTyphi scheme, we needed to create fake levels within the scheme to facilitate correct calling by Mykrobe. These fake levels are not reported in the output generated by the parser, but they are present in the raw json files output by Mykrobe. These can be identified in the json output as they are prepended by the word "lineage", and will always have a call of 0 from Mykrobe.

**Dependencies**

* Python 3 
* library `pandas`

**Input**
* json files output by Mykrobe
* prefix for output file

**Output**
* tsv file, one row per json file

**Example command**
```
python parse_typhi_mykrobe.py --jsons *.json --prefix mykrobe_out
```

### Explanation of columns in the output:
* **genome**: sample ID
* **species**: species call from Mykrobe (typhi or unknown; determined by matching to Typhi STs from the 7-locus MLST scheme)
* **spp_percent:** percentage coverage to the Typhi ST probes
* **final genotype**: final genotype call from Mykrobe, using the scheme currently implemented in Genotyphi.
* **confidence**: measure of confidence in the final genotype call, summarising read support across all levels in the hierarchy (lineage, clade, subclade, etc)
  * _strong_ - high quality calls (quality of '1' reported by Mykrobe) for ALL levels;
  * _moderate_ - reduced confidence for ONE node (Mykrobe quality of '0.5', but still with >50% read support for the marker allele), high quality calls for ALL OTHERS;
  * _weak_ - low quality for one or more nodes (Mykrobe quality of '0.5' and <50% read support OR Mykrobe quality of '0').
* **acrB_R717L, acrB_7171Q**: calls for individual acrB (azithromycin resistance) mutation. 0 indicates mutation is absent, 1 indicates mutation is present.
* **num QRDR**: Total number of mutations detected in the quinolone-resistance determining regions (QRDR) of genes _gyrA_, _parC_ and _gyrB_.
* **lowest support for genotype marker**: For any markers in the final genotype call that do not have a Mykrobe quality of '1', this column reports the percentage of reads supporting the marker allele at the most poorly supported marker (details of all such markers appear in the 'poorly supported markers column').
* **poorly supported markers**: Lists any markers in the final genotype call that do not have Mykrobe quality of '1'. Markers are separated by ';', values in brackets represent the quality call from Mykrobe, followed by the read depth at the alternate / reference alleles. The lowest read support amongst these markers is reported in the previous column.
* **max support for additional markers**: For any markers detected that are incongruent with the final genotype call, this column reports the percentage of reads supporting the marker allele at the best supported additional marker.
* **additional markers**: Lists any markers that are incongruent with the final genotype call. Markers are separated by ';', and the format is identical to column _poorly supported markers_. The highest read support for any such marker is reported in the previous column.
* **node support**: A list of all markers in the final genotype call with their Mykrobe quality calls (1, 0.5, or 0) and the read depths at the marker allele / reference allele (as per _poorly supported markers_ and _additional markers_).
* **parC_S80R, parC_S80I, parC_E84G, parC_E84K, gyrA_S83F, gyrA_S83Y, gyrA_D87G, gyrA_D87N, gyrA_D87V, gyrA_D87Y, gyrB_S464F, gyrB_S464Y**: calls for each individual QRDR mutation. 0 indicates mutation is absent, 1 indicates mutation is present.
* **remaining columns** indicate presence (1) or absence (0) of each AMR gene or plasmid replicon as indicated by the header of the column. For the AMR genes, _specific alleles are **not detected**_ (with the exception of dfrA, where Mykrobe can distinguish between dfrA7, dfrA5 and dfrA15).
  * **IncHI1_ST6** indicates whether the IncHI1 replicon (from column IncHI1BR27) is pST6 (1) or not (0).


## Original Implementation: Data must be pre-mapped to CT18 reference, inputs are BAMs or VCFs

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
