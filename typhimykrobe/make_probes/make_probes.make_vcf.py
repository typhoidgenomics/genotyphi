#!/usr/bin/env python3

import glob
import logging
from cluster_vcf_records import vcf_file_read, vcf_record

log = logging.getLogger()
log.setLevel(logging.INFO)

def load_variants_from_vcf(filename):
    logging.info(f"load VCF {filename}")
    header, records = vcf_file_read.vcf_file_to_list(filename)
    variants = set()
    for record in records:
        # ignore anything that isn't a SNP
        if len(record.REF) != 1 or any([len(x) != 1 for x in record.ALT]):
            continue

        # Only use records where we have the genotype and it's non-ref
        # homozygous
        if "GT" not in record.FORMAT:
            continue

        gt_set = set(record.FORMAT["GT"].split("/"))
        if "." in gt_set or len(gt_set) != 1 or "0" in gt_set:
            continue

        gt = int(gt_set.pop())
        alt = record.ALT[gt-1]
        variants.add((record.CHROM, record.POS, record.REF, alt))

    return variants


def load_all_variants(vcf_files):
    variants = set()
    for vcf_file in vcf_files:
        variants.update(load_variants_from_vcf(vcf_file))
    variants = list(variants)
    variants.sort()
    return variants


def write_vcf(variants, outfile):
    with open(outfile, "w") as f:
        print("##fileformat=VCFv4.2", file=f)
        print('##FORMAT=<ID=GT_CONF,Number=1,Type=Float,Description="Genotype confidence. Difference in log likelihood of most likely and next most likely genotype">', file=f)
        print("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "sample_name", sep="\t", file=f)
        for i, (chrom, pos, ref, alt) in enumerate(variants):
            print(chrom, pos+1, i, ref, alt, "255", "PASS", "SVTYPE=SNP", "GT:GT_CONF", "1/1:100", sep="\t", file=f)


vcf_files = glob.glob("vcfs_for_probe_creation/*.vcf")
variants = load_all_variants(vcf_files)
logging.info("Loaded all variants. Writing VCF file")
write_vcf(variants, "make_probes.backgrounds.vcf")
logging.info("Finished writing VCF file")


