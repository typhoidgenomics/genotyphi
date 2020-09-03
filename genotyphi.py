#!/usr/bin/env python

"""
Input BAM (recommended) or VCF (if highly trusted SNP data) relative to Typhi CT18 (AL513382) and assign Typhi genotype codes and detect AMR mutations.

Authors - Kat Holt (kholt@unimelb.edu.au), Zoe Dyson (zoe.dyson@lshtm.ac.uk & zad24@medschl.cam.ac.uk)

Documentation - https://github.com/katholt/genotyphi

Dependencies:
 SAMtools (v1.2) and bcftools (v1.2) are required to genotype from BAMs.

Last modified - Sep 3rd, 2020
Adapted from https://github.com/katholt/genotyphi/commit/55eb40b71390ee5a39f2cc51fd95796f5ee26e65
"""

import argparse
import gzip
import logging
import os
import pathlib
import re
import sys
import tempfile
from subprocess import call, Popen, PIPE
from typing import List, Optional

from Bio.SeqIO.FastaIO import SimpleFastaParser

LOG_FORMAT = '\033[2;36m{asctime}\033[0;1m {levelname:<8}\033[0m {message} \033[2m[{name}:{funcName}:{lineno}]\033[0m'


def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description='VCF to Typhi genotypes')
    parser.add_argument('--mode', choices=('vcf', 'bam', 'vcf_parsnp'), default='bam',
                        help='Mode to run in based on input files (vcf, bam, or vcf_parsnp)')
    parser.add_argument('--vcf', nargs='+', type=str, required=False,
                        help='VCF file(s) to genotype (Mapping MUST have been done using CT18 as a reference sequence)')
    parser.add_argument('--bam', nargs='+', type=str, required=False,
                        help='BAM file(s) to genotype (Mapping MUST have been done using CT18 as a reference sequence)')
    parser.add_argument('--ref_id', type=str, required=False,
                        help='Name of the reference in the VCF file (#CHROM column) or fasta file. Note that CT18 has '
                             'genotype 3.2.1. If all your strains return this genotype, it is likely you have '
                             'specified the name of the refrence sequence incorrectly; please check your VCFs.')
    parser.add_argument('--phred', type=int, required=False, default=20,
                        help='Minimum phred quality to count a variant call vs CT18 as a true SNP (default 20)')
    parser.add_argument('--min_prop', type=float, required=False, default=0.1,
                        help='Minimum proportion of reads required to call a SNP (default 0.1)')
    parser.add_argument('--ref', type=str, required=False,
                        help='Reference sequence in fasta format. Required if bam files provided.')
    parser.add_argument('--output', type=str, required=False, default=None,
                        help='Location and name for output file (default=stdout)')
    return parser.parse_args()


### Genotype SNP definitions

loci = [655112, 773487, 1804415, 1840727, 3640678, 270120, 102135, 316489, 4105384, 555826, 2360997, 4664137, 2166082,
        30192, 4288272, 2737027, 1215983, 4132985, 146673, 2517324, 3920009, 3276735, 1173984, 2683312, 4013386,
        4094437, 827036, 3476114, 4355243, 2723847, 4388609, 703762, 431216, 3095443, 316186, 1934711, 2811222, 3092900,
        2723724, 3437570, 1780319, 1535365, 1792810, 3062270, 1799842, 432732, 3069182, 2732615, 3770391, 2269835,
        4215341, 4602946, 3368641, 2245432, 3164162, 3923165, 1811809, 3729635, 3817752, 183033, 1615350, 2342045,
        3996717, 2640029, 989024, 3806278, 1611156, 2348902, 1193220, 3694947, 955875,
        3498544,
        2424394,
        2272144,
        561056, 3164873]
snp_alleles = ['T', 'A', 'C', 'A', 'A', 'T', 'A', 'C', 'T', 'A', 'T', 'T', 'A', 'T', 'T', 'G', 'G', 'A', 'A', 'A', 'T',
               'A', 'T', 'A', 'A', 'A', 'A', 'T', 'C', 'A', 'C', 'A', 'A', 'A', 'T', 'T', 'T', 'T', 'G', 'C', 'A', 'T',
               'T', 'C', 'C', 'T', 'G', 'A', 'T', 'G', 'C', 'T', 'T', 'A', 'A', 'A', 'T', 'T', 'T', 'T', 'T', 'A', 'T',
               'A', 'T', 'A', 'A', 'T', 'C', 'G', 'A',
               'G',
               'A',
               'A',
               'A', 'A']
groups = ['0.1', '0.0.1', '0.0.2', '0.0.3', '0.1.1', '0.1.2', '0.1.3', '1', '1.1', '1.1.1', '1.1.2', '1.1.3', '1.1.4',
          '1.2', '1.2.1', '2', '2.0.1', '2.0.2', '2.1', '2.1.1', '2.1.2', '2.1.3', '2.1.4', '2.1.5', '2.1.6', '2.1.7',
          '2.1.8', '2.1.9', '2.2', '2.2.1', '2.2.2', '2.2.3', '2.2.4', '2.3', '2.3.1', '2.3.2', '2.3.3', '2.3.4',
          '2.3.5', '2.4', '2.4.1', '2.5', '2.5.1', '3', '3.0.1', '3.0.2', '3.1', '3.1.1', '3.1.2', '3.2', '3.2.1',
          '3.2.2', '3.3', '3.3.1', '3.4', '3.5', '3.5.1', '3.5.2', '3.5.3', '3.5.4', '4', '4.1', '4.1.1', '4.2',
          '4.2.1', '4.2.2', '4.2.3', '4.3.1', '4.3.1.1', '4.3.1.2', '4.3.1.1.P1',
          '3.3.2',
          '3.3.2.Bd1',
          '3.3.2.Bd2',
          '4.3.1.3', '2.5.2']

### QRDR SNP definitions

qrdr_loci = [523109, 523109, 2333762, 2333762, 2333750, 2333750, 2333751, 2333751, 3196469, 3196470, 3196458, 3196459]
qrdr_snp_alleles = ['A', 'T', 'A', 'T', 'A', 'C', 'A', 'T', 'T', 'A', 'C', 'T']
qrdr_groups = [' acrB-R717L', ' acrB-R717Q', ' gyrA-S83F', ' gyrA-S83Y', ' gyrA-D87V', ' gyrA-D87G', ' gyrA-D87Y',
               ' gyrA-D87N', ' parC-S80R', ' parC-S80I', ' parC-E84G', ' parC-E84K']


# check if this SNP defines a QRDR group
def checkQRDRSNP(vcf_line_split, this_qrdr_groups, args):
    qrdr_snp = int(vcf_line_split[1])
    if qrdr_snp in qrdr_loci and float(vcf_line_split[5]) > args.phred:
        logging.info(vcf_line_split)
        m = re.search(r"DP4=(\d+),(\d+),(\d+),(\d+)", vcf_line_split[7])
        if m is None:
            if vcf_line_split[4] == '.':
                qrdr_snp_proportion = float(
                    -1)  # set unknowns to negative so that we know this is not a real proportion

            else:  # if the ALT is not '.' i.e. if the alt is not same as ref
                try:
                    ad = vcf_line_split[9].split(':')[1].split(',')  # get the AD ratio
                    alt_read_count = int(ad[1])
                    total_read_count = int(ad[0]) + alt_read_count
                    qrdr_snp_proportion = float(alt_read_count) / total_read_count
                except IndexError:
                    qrdr_snp_proportion = float(-1)

        else:
            alt_read_count = int(m.group(3)) + int(m.group(4))
            total_read_count = alt_read_count + int(m.group(1)) + int(m.group(2))
            if float(total_read_count) == 0:
                qrdr_snp_proportion = float(-1)
            else:
                qrdr_snp_proportion = float(alt_read_count) / total_read_count
        if qrdr_snp_proportion > args.min_prop:
            qrdr_snp_allele = vcf_line_split[4]
            for idx, locus in enumerate(qrdr_loci):
                if qrdr_snp == locus and qrdr_snp_allele == qrdr_snp_alleles[idx]:
                    this_qrdr_groups.append(qrdr_groups[idx])  # Add QRDR SNP

    return this_qrdr_groups


# check if this SNP defines a group
def checkSNP(vcf_line_split, this_groups, proportions, args):
    location = int(vcf_line_split[1])
    if location in loci:
        i = loci.index(location)

        if float(vcf_line_split[5]) > args.phred:
            logging.debug(vcf_line_split)
            m = re.search(r"DP4=(\d+),(\d+),(\d+),(\d+)", vcf_line_split[7])
            if m is None:
                if vcf_line_split[4] == '.':
                    snp_proportion = float(-1)  # set unknowns to negative so that we know this is not a real proportion
                else:  # if the ALT is not '.' i.e. if the alt is not same as ref
                    try:
                        ad = vcf_line_split[9].split(':')[1].split(',')  # get the AD ratio
                        alt_read_count = int(ad[1])
                        total_read_count = int(ad[0]) + alt_read_count
                        snp_proportion = float(alt_read_count) / total_read_count
                    except IndexError:
                        snp_proportion = float(-1)

            else:
                alt_read_count = int(m.group(3)) + int(m.group(4))
                total_read_count = alt_read_count + int(m.group(1)) + int(m.group(2))
                if float(total_read_count) == 0:
                    snp_proportion = float(-1)
                else:
                    snp_proportion = float(alt_read_count) / total_read_count
            if snp_proportion > args.min_prop:
                this_allele = vcf_line_split[4]
                if this_allele == snp_alleles[i]:
                    this_groups.append(groups[i])  # retrieve the group that this SNP defines
                    proportions[groups[i]] = snp_proportion
    return this_groups, proportions


def checkSNPmulti(vcf_line_split, this_groups, args):
    snp = int(vcf_line_split[1])
    if snp in loci:
        i = loci.index(snp)
        strain = 0
        for gt in vcf_line_split[10:]:
            if (int(gt) == 1) and (vcf_line_split[4] == snp_alleles[i]):
                if strain in this_groups:
                    this_groups[strain].append(groups[i])  # retrieve the group that this SNP defines
                else:
                    this_groups[strain] = [groups[i]]
            strain += 1
    return this_groups


def parseGeno(this_groups, proportions) -> str:
    """sort groups into the three levels (primary, clade, subclade)"""
    subclades = []
    clades = []
    primary = []
    for group in this_groups:
        level = len(group.split("."))
        if level == 5:
            subclades.append(group)
        if level == 4:
            subclades.append(group)
        if level == 3:
            subclades.append(group)
        elif level == 2:
            clades.append(group)
        elif level == 1:
            primary.append(group)

    # fix 4.3.1/4.3.1.1/4.3.1.2/4.3.1.P1/4.3.1.3 nesting
    if ('4.3.1.3' in subclades) and ('4.3.1' in subclades):
        subclades.remove('4.3.1')
    if ('4.3.1.1' in subclades) and ('4.3.1' in subclades):
        subclades.remove('4.3.1')
    if ('4.3.1.2' in subclades) and ('4.3.1' in subclades):
        subclades.remove('4.3.1')
    if ('4.3.1.1.P1' in subclades) and ('4.3.1' in subclades):
        subclades.remove('4.3.1')
    if ('4.3.1.1.P1' in subclades) and ('4.3.1.1' in subclades):
        subclades.remove('4.3.1.1')

    # fix 3.3.2.Bd nesting
    if ('3.3.2.Bd1' in subclades) and ('3.3.2' in subclades):
        subclades.remove('3.3.2')
    if ('3.3.2.Bd2' in subclades) and ('3.3.2' in subclades):
        subclades.remove('3.3.2')

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
            primary = ['2']  # clade 2 differs from CT18 by the clade3-defining SNP
        # note other option is clade 4 snp, which defines clade 4 relative to CT18
    elif len(primary) == 2:
        if ('2' in primary) and ('3' in primary):
            primary = ['1']  # clade 2 differs from CT18 by the clade3-defining SNP
    elif len(primary) == 0:
        primary = ['3']
    elif len(primary) == 3:
        if ('1' in primary) and ('2' in primary) and ('3' in primary):
            primary = ['0']

    # fix clade relative to CT18:
    if '3.2' in clades:
        clades.remove('3.2')  # anything NOT in 3.2 will have this SNP
    else:
        if len(clades) == 0:
            clades.append('3.2')  # anything with no clade, and 3.2 SNP not called, belongs in 3.2 with CT18

    # fix 3.5.3, 3.5.4 nesting
    if ('3.5.3' in clades) and ('3.5.4' in clades):
        clades.remove('3.5.3')

    # fix subclades relative to CT18:
    if '3.2.1' in subclades:
        subclades.remove('3.2.1')  # anything NOT in 3.2.1 will have this SNP
    else:
        if len(subclades) == 0:
            subclades.append('3.2.1')  # anything with no subclade, and 3.2.1 SNP NOT called, belongs in 3.2.1 with CT18

        # add zero-th clade/subclade where unresolved -- disabled
        # if len(clades) == 0:
        #	if len(primary) == 1:
        #		clades.append(primary[0] + '.0')
        # if len(subclades) == 0:
    # if len(clades) == 1:
    #		subclades.append(clades[0] + '.0')

    # store final genotype, to the lowest level available
    final_geno = primary[0]
    if len(clades) > 0:
        final_geno = ','.join(clades)
    if len(subclades) > 0:
        final_geno = ','.join(subclades)

    # add proportion of reads supporting each of these groups
    p_prod = 1

    p_sub = []
    for group in subclades:
        if group in proportions:
            p_sub.append(str(round(proportions[group], 2)))
            p_prod = p_prod * proportions[group]

    p_cl = []
    for group in clades:
        if group in proportions:
            p_cl.append(str(round(proportions[group], 2)))
            p_prod = p_prod * proportions[group]

    p_pr = []
    for group in primary:
        if group in proportions:
            p_pr.append(str(round(proportions[group], 2)))
            p_prod = p_prod * proportions[group]

    # final call
    info = final_geno + '\t'
    if 'A' in proportions:
        info += 'A'  # annotate as 'A' to indicate this comes from assembled data and not reads
    else:
        info += str(round(p_prod, 2))  # indicate proportion of reads supporting this call

    # level calls
    info += '\t' + ','.join(subclades) + '\t' + ','.join(clades) + '\t' + ','.join(primary)

    # level proportions
    info += '\t' + ','.join(p_sub) + '\t' + ','.join(p_cl) + '\t' + ','.join(p_pr)

    return info


# exception to raise if the command we try to run fails for some reason
class CommandError(Exception):
    pass


def run_command(command, **kwargs):
    """Execute a shell command and check the exit status and any O/S exceptions"""
    command_str = ' '.join(command)
    logging.info('Running: {}'.format(command_str))
    try:
        exit_status = call(command, **kwargs)
    except OSError as e:
        message = "Command '{}' failed due to O/S error: {}".format(command_str, str(e))
        raise CommandError({"message": message})
    if exit_status != 0:
        message = "Command '{}' failed with non-zero exit status: {}".format(command_str, exit_status)
        raise CommandError({"message": message})


# main function
def main():
    args = parse_args()

    logging.basicConfig(format=LOG_FORMAT, level=logging.DEBUG, style='{')
    logging.info(f'Parsed args: {args}')
    ref_id = args.ref_id
    ref_fasta = args.ref
    bams: Optional[List[str]] = args.bam
    vcfs: Optional[List[str]] = args.vcf
    mode: str = args.mode

    if (((mode == 'vcf') and vcfs and ref_id) or (
            (mode == 'bam') and bams and ref_fasta and ref_id) or (
            (mode == 'vcf_parsnp') and vcfs)):

        if mode == 'bam':
            # GENERATE VCFS (1 per strain) FROM BAMS
            vcfs = vcfs_from_bams(bams, ref_fasta, ref_id, phred=args.phred)

        output_headers = get_output_headers(mode)
        # PARSE MAPPING BASED VCFS (1 per strain)
        if vcfs and (mode != 'vcf_parsnp'):
            with open(args.output, 'w') if args.output else sys.stdout as output_file:
                output_file.write(f'{output_headers}\n')
                for vcf in vcfs:
                    snp_count = 0
                    this_groups = []  # list of groups identified by defining SNPs
                    this_qrdr_groups = []  # list of QRDR SNPs found
                    # proportion of reads supporting each defining SNP; key = group, value = proportion
                    proportions = {}
                    _, ext = os.path.splitext(vcf)
                    with gzip.open(vcf, 'r') if ext == '.gz' else open(vcf, 'r') as vcf_handle:
                        any_ref_line = 0
                        for line in vcf_handle:
                            if not line.startswith('#'):
                                vcf_line_split = line.strip().split('\t')
                                if not vcf_line_split:
                                    continue
                                snp_count += 1
                                if vcf_line_split[0] == args.ref_id:
                                    # parse this SNP line
                                    any_ref_line = 1
                                    (this_groups, proportions) = checkSNP(vcf_line_split, this_groups, proportions,
                                                                          args)
                                    this_qrdr_groups = checkQRDRSNP(vcf_line_split, this_qrdr_groups, args)

                    logging.info("qrdr groups".join(qrdr_groups))
                    if any_ref_line > 0:
                        genotyping_result: str = parseGeno(this_groups, proportions)
                        if args.bam:
                            output_file.write(f'{vcf}'
                                              f'\t{genotyping_result}'
                                              f'\t{snp_count}'
                                              f'\t{",".join(this_qrdr_groups)}'
                                              f'\n')
                        else:
                            output_file.write(f'{vcf}'
                                              f'\t{genotyping_result}'
                                              f'\t{",".join(this_qrdr_groups)}'
                                              f'\n')
                    else:
                        output_file.write(f'{vcf}'
                                          f'\tNo SNPs encountered against expected reference. '
                                          f'Wrong reference or no SNP calls?'
                                          f'\n')

        # PARSE PARSNP VCF (multiple strains)

        if mode == 'vcf_parsnp':

            if not args.ref_id:
                args.ref_id = '1'

            for vcfm in args.vcf:

                # read file
                (file_name, ext) = os.path.splitext(vcfm)

                if ext == '.gz':
                    vcf_handle = gzip.open(vcfm, 'r')
                else:
                    vcf_handle = open(vcfm, 'r')

                any_ref_line = 0

                this_groups = {}  # list of groups identified by defining SNPs, key = strain id (column number)
                strains = []

                for line in vcf_handle:
                    vcf_line_split = line.rstrip().split()
                    if vcf_line_split[0] == '#CHROM':
                        strains = vcf_line_split[10:]
                    if not line.startswith('#'):
                        if vcf_line_split[0] == args.ref_id:
                            any_ref_line = 1  # parse this SNP line
                            this_groups = checkSNPmulti(vcf_line_split, this_groups, args)

                vcf_handle.close()

                # collate by strain
                if any_ref_line > 0:
                    for strain in this_groups:
                        genotyping_result = parseGeno(this_groups[strain], ['A'])
                        output_file.write(strains[strain] + '\t' + genotyping_result + '\n')
                else:
                    output_file.write(f'{strains[strain]}'
                                      f'\tNo SNPs encountered against expected reference. '
                                      f'Wrong reference or no SNP calls?\n')

        output_file.close()
    else:
        logging.info('Missing or incomplete input parameters, please check these and try again.')


def get_output_headers(mode):
    output_header = ['File', 'Final_call', 'Final_call_support', 'Subclade', 'Clade', 'PrimaryClade',
                     'Support_Subclade', 'Support_Clade', 'Support_PrimaryClade']
    if mode == 'bam':
        output_header += ['Number of SNPs called', 'AMR mutations']
    elif mode == 'vcf':
        output_header += ['AMR mutations']
    output_header = '\t'.join(output_header)
    return output_header


def vcfs_from_bams(bams, ref_fasta, ref_id, phred=20):
    vcfs = []
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = pathlib.Path(tmpdir)
        temp_fasta_path = tmpdir_path / 'tmp.fasta'
        create_faidx_fasta(ref_fasta, ref_id, temp_fasta_path)
        # index fasta file
        run_command(['samtools', 'faidx', str(temp_fasta_path)])
        fai_path = pathlib.Path(f'{temp_fasta_path}.fai')
        assert fai_path.exists(), f'Reference FASTA index file "{fai_path}" does not exist.'
        bed_path = tmpdir_path / 'tmp.bed'
        create_bed_file(ref_id, bed_path)
        create_bed_file(ref_id, f'genotyphi.bed')
        for bam in bams:
            bam = pathlib.Path(bam)
            vcf_path = pathlib.Path(f'{bam}.vcf')
            logging.info(f'bam files supplied, generating vcf file for {bam}')
            if not pathlib.Path(f'{bam}.bai').exists():  # index bam file if indexed bam not provided
                run_command(['samtools', 'index', str(bam)])
            samtools_mpileup_cmd = ['samtools', 'mpileup', '-q', str(phred), '-ugB', '-f', str(temp_fasta_path), '-l',
                                    str(bed_path), '-I', str(bam)]
            logging.info(f'Running {" ".join(samtools_mpileup_cmd)}')
            proc_samtools_mpileup = Popen(samtools_mpileup_cmd, stdout=PIPE)
            p = Popen(['bcftools', 'call', '-c', '-o', str(vcf_path)], stdin=proc_samtools_mpileup.stdout)
            p.communicate()
            assert vcf_path.exists(), f'VCF file "{vcf_path}" does not exist for BAM "{bam}"'
            vcfs.append(vcf_path)
    return vcfs


def create_bed_file(ref_id, bed_path) -> None:
    # coordinates in zero-base, half-open for SAMtools compatible bed file
    # create temporary bed file for SAMtools
    with open(bed_path, 'w') as temp_bed_file:
        for locus in sorted(loci + qrdr_loci):
            temp_bed_file.write(f'{ref_id}\t{locus - 1}\t{locus}\n')


def create_faidx_fasta(ref_fasta, ref_id, temp_fasta_path) -> None:
    # create SAMtools compatible fasta file
    with open(ref_fasta, 'r') as fasta_file, open(temp_fasta_path, 'w') as temp_fasta_file:
        for header, seq in SimpleFastaParser(fasta_file):
            if ref_id in header:
                temp_fasta_file.write(f'>{ref_id}\n{seq}\n')


# call main function
if __name__ == '__main__':
    main()
