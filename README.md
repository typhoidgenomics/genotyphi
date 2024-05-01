# Genotyping Salmonella Typhi

This repository houses the GenoTyphi genotyping scheme for Salmonella Typhi. It also describes how to call genotypes, AMR and plasmid markers in Typhi sequence read sets using Mykrobe ('Typhi Mykrobe') and links to alternative tools for calling genotypes from reads or assemblies.

* [GenoTyphi scheme](https://github.com/typhoidgenomics/genotyphi#GenoTyphi_Typing_Scheme)
* [Typing from reads using Mykrobe ('Typhi Mykrobe')](https://github.com/typhoidgenomics/genotyphi#Typhi_Mykrobe)
* [Other tools for callling genotypes from reads or assemblies](https://github.com/typhoidgenomics/genotyphi#Other_typing_tools)


## GenoTyphi Scheme

The GenoTyphi genotyping scheme divides the *Salmonella* Typhi population into 4 major lineages, and >75 different clades and subclades. The scheme specification is detailed in the file `GenoTyphi_specification.csv` in this repository, which also includes the standard clade-level colour codes that we use for consistency across papers.

The scheme is described in this paper, **which serves as the primary citation for the scheme**: ["An extended genotyping framework for Salmonella enterica serovar Typhi, the cause of human typhoid", Wong et al, 2016, Nature Communications](http://www.nature.com/articles/ncomms12827/).

Subsequent updates to the genotyping scheme, including new genotypes and mutations conferring resistance to fluoroquinolones and azithromycin, are summarised in ["Five years of GenoTyphi: updates to the global Salmonella Typhi genotyping framework", Dyson & Holt, 2021, Journal of Infectious Diseases](https://doi.org/10.1093/infdis/jiab414).

Two implementations of the code that can assign genotypes to Typhi genomes are available in this repository:

* [Mykrobe implementation](https://github.com/katholt/genotyphi/blob/main/README.md#mykrobe-implementation), which takes as input **sequence reads (fastq files)** from Illumina or long-read platforms. It also detects known mutations in the quinolone-resistance determining region (QRDR) of genes *gyrA*, *gyrB* and *parC*; and the *acrB*-R717Q/L mutations associated with azithromycin resistance; acquired antimicrobial resistance genes; plasmid replicons and major subtypes of the IncHI1 plasmid typically associated with multidrug resistance (list of targets is in 'AMR_genes_mutations_plasmids.csv'). Drugs for which resistance is typed are: `ampicillin`, `azithromycin`, `ceftriaxone`, `ciprofloxacin`, `chloramphenicol`, `sulfonamides`, `trimethoprim`, `trimethoprim-sulfamethoxazole`, `tetracycline` indicate resistant (R) or susceptible (S) predictions for each genome.

**Whichever tool you use to access the GenoTyphi scheme, please cite the [GenoTyphi paper](https://doi.org/10.1093/infdis/jiab414).**

**If you use the scripts in this repository, please also cite the repository:** [![DOI](https://zenodo.org/badge/45819844.svg)](https://zenodo.org/badge/latestdoi/45819844)

## Typhi Mykrobe

### Quick start

#### Install Mykrobe:

From bioconda:
```
conda install -c bioconda mykrobe
```
From source (after downloading mykrobe, from the mykrobe directory):
```
pip3 install . && mykrobe panels update_metadata && mykrobe panels update_species all
```

#### Run Mykrobe on fastq file/s for a given genome:

```
mykrobe predict --sample aSample \
  --species typhi \
  --format json \
  --out aSample.json \
  --seq aSample_1.fastq.gz aSample_2.fastq.gz
```

Output is one JSON file per genome

#### Tabulate Mykrobe results for one or more genomes:

(requires Python3 + pandas library)

```
python parse_typhi_mykrobe.py --jsons *.json --prefix mykrobe_out
```

Output is a single tab-delimited table, output format is [described below](#explanation-of-columns-in-the-output).

### Detailed instructions
* [Install Mykrobe and Typhi panels](#installing-mykrobe)
* [Run Mykrobe](#running-mykrobe)
* [Tabulate Mykrobe output](#parse-mykrobe-output)
* [Output format](#explanation-of-columns-in-the-output)

### Installing mykrobe

First, install Mykrobe (v0.10.0+) as per the instructions on the [Mykrobe github](https://github.com/Mykrobe-tools/mykrobe).

Once Mykrobe is installed, you can run the following two commands to ensure you have the most up-to-date panels for genotyping, including the [Typhi panel](https://doi.org/10.6084/m9.figshare.21695528.v1):
```
mykrobe panels update_metadata
mykrobe panels update_species all
```

You can check what version of the scheme is currently loaded in your Mykrobe installation via:
```
mykrobe panels describe
```

### Running Mykrobe

Inputs are fastq files.

Mykrobe can be run on each individual sample using the command below. Replace `aSample` with the name of your isolate. The command below uses Illumina data.

```
mykrobe predict --sample aSample --species typhi --format json --out aSample.json --seq aSample_1.fastq.gz aSample_2.fastq.gz
```

To run on ONT data instead, add the `--ont` flag to your command.

Further details on options can be found on the Mykrobe wiki: https://github.com/Mykrobe-tools/mykrobe/wiki

### Parse Mykrobe output

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
* **serovar**: species call from Mykrobe (Typhi or unknown; determined by matching to Typhi STs from the 7-locus MLST scheme)
* **spp_percent:** percentage coverage to the Typhi ST probes
* **final genotype**: final genotype call from Mykrobe, using the latest version of GenoTyphi
* **confidence**: measure of confidence in the final genotype call, summarising read support across all levels in the hierarchy (lineage, clade, subclade, etc)
  * _strong_ - high quality calls (quality of '1' reported by Mykrobe) for ALL levels;
  * _moderate_ - reduced confidence for ONE node (Mykrobe quality of '0.5', but still with >50% read support for the marker allele), high quality calls for ALL OTHERS;
  * _weak_ - low quality for one or more nodes (Mykrobe quality of '0.5' and <50% read support OR Mykrobe quality of '0').
* **lowest support for genotype marker**: For any markers in the final genotype call that do not have a Mykrobe quality of '1', this column reports the percentage of reads supporting the marker allele at the most poorly supported marker (details of all such markers appear in the 'poorly supported markers column').
* **poorly supported markers**: Lists any markers in the final genotype call that do not have Mykrobe quality of '1'. Markers are separated by ';', values in brackets represent the quality call from Mykrobe, followed by the read depth at the alternate / reference alleles. The lowest read support amongst these markers is reported in the previous column.
* **max support for additional markers**: For any markers detected that are incongruent with the final genotype call, this column reports the percentage of reads supporting the marker allele at the best supported additional marker.
* **additional markers**: Lists any markers that are incongruent with the final genotype call. Markers are separated by ';', and the format is identical to column _poorly supported markers_. The highest read support for any such marker is reported in the previous column.
* **node support**: A list of all markers in the final genotype call with their Mykrobe quality calls (1, 0.5, or 0) and the read depths at the marker allele / reference allele (as per _poorly supported markers_ and _additional markers_).
* **resistance predictions**: Columns `ampicillin`, `azithromycin`, `ceftriaxone`, `ciprofloxacin`, `chloramphenicol`, `sulfonamides`, `trimethoprim`, `trimethoprim-sulfamethoxazole`, `tetracycline` indicate resistant (R) or susceptible (S) predictions for each genome. A resistant prediction is appended with the genes or mutations giving rise to that prediction.
* **remaining columns** indicate presence (1) or absence (0) of each QRDR or _acrB_ mutation, AMR gene, or plasmid replicon as indicated by the header of the column. For the AMR genes, _specific alleles are **not detected**_ (with the exception of _dfrA_ and blaOXA, where Mykrobe can distinguish between _dfrA1_, _dfrA5_, _dfrA7_, _dfrA14_, _dfrA15_, _dfrA17_ and _dfrA18_; and blaOXA-7 and blaOXA-134).
  * **num QRDR**: Total number of mutations detected in the quinolone-resistance determining regions (QRDR) of genes _gyrA_, _parC_ and _gyrB_.
  * **IncHI1_ST6** indicates whether the IncHI1 replicon (from column IncHI1BR27) is pST6 (1) or not (0).
 
Sequences and details of probes are available [here](https://doi.org/10.26180/14881701.v1).


## Other typing tools
* [Original Python implementation](https://github.com/katholt/genotyphi/blob/main/README.md#original-implementation-pre-mapped-data), which takes as input **BAM or VCF files that the user has already generated** by mapping Illumina reads to the reference genome CT18. It also detects the QRDR and *acrB* mutations listed in 'AMR_genes_mutations_plasmids.csv' in this repository.

The GenoTyphi scheme is also available via the online analysis platform [Pathogenwatch](https://pathogen.watch/), which facilitates automated analysis of Typhi genome assemblies as described in [this paper](https://www.nature.com/articles/s41467-021-23091-2).

