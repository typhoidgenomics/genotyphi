#!/usr/bin/env python3

import argparse
import csv
import json
import os
import pyfastaq


def load_profile_file(filename, wanted_sts):
    wanted_seqs = {}
    wanted_profiles = {}
    with open(filename) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for d in reader:
            if d["ST"] not in wanted_sts:
                continue

            st = int(d["ST"])
            del d["ST"]
            #del d["clonal_complex"]
            wanted_profiles[st] = {k: int(v) for k, v in d.items()}
            for k, v in d.items():
                if k not in wanted_seqs:
                    wanted_seqs[k] = set()
                wanted_seqs[k].add(int(v))
    return wanted_seqs, wanted_profiles


def load_wanted_seqs(wanted_seqs, fasta_dir):
    seqs = {}
    for gene, alleles in wanted_seqs.items():
        seqs[gene] = {}
        reader = pyfastaq.sequences.file_reader(os.path.join(fasta_dir, f"{gene}.fas"))
        for seq in reader:
            seq_gene, seq_number = seq.id.split("_")
            seq_number = int(seq_number)
            assert seq_gene == gene
            if seq_number in alleles:
                seqs[gene][seq_number] = seq.seq

    return seqs


def check_got_all_seqs(wanted_seqs, got_seqs):
    for gene, alleles in wanted_seqs.items():
        assert gene in got_seqs
        for allele in alleles:
            assert allele in got_seqs[gene]


def write_probes_fasta(wanted_profiles, seqs, outfile, seq_name):
    data = {}
    with open(outfile, "w") as f:
        for st, profile in wanted_profiles.items():
            seqs_to_print = []
            for gene, allele in profile.items():
                seqs_to_print.append(seqs[gene][allele])
            sequence = "N".join(seqs_to_print)
            name = "&".join([
                f"st_{st}?name={seq_name}",
                f"length={len(sequence)}",
                "panel_type=species",
            ])
            seq = pyfastaq.sequences.Fasta(name, sequence)
            data[st] = {"probe_name": name, "profile": profile}
            print(seq, file=f)

    return data


parser = argparse.ArgumentParser(
    description = "Given ariba pubmlst_download directory, makes mykrobe probes for given clonal complex",
    usage="%(prog)s <sts> <mlst_dir> <seq_name> <outprefix>",
)

parser.add_argument("st_list", help="File listing the STs to include")
parser.add_argument("mlst_dir", help="Name of input directory with the allele sequences and profile table from pubMLST")
parser.add_argument("seq_name", help="Name in the 'name=...' part of fasta header lines")
parser.add_argument("outprefix", help="Prefix of output files")
options = parser.parse_args()

#fasta_dir = os.path.join(options.ariba_dir, "pubmlst_download")
fasta_dir = options.mlst_dir
profile_txt = os.path.join(fasta_dir, "profile.txt")
# get a list of all the STs we want to include from the file
wanted_sts = []
with open(options.st_list, 'r') as f:
    for line in f:
        wanted_sts.append(line.strip())
# pull out dictionaries of the wanted sequences and profiles for the STs of interest
wanted_seqs, wanted_profiles = load_profile_file(profile_txt, wanted_sts)
#print('wanted seqs')
#print(wanted_seqs)
#print('wanted profiles')
#print(wanted_profiles)
seqs = load_wanted_seqs(wanted_seqs, fasta_dir)
check_got_all_seqs(wanted_seqs, seqs)
data = write_probes_fasta(wanted_profiles, seqs, options.outprefix + ".fasta", options.seq_name)
with open(options.outprefix + ".json", "w") as f:
    json.dump(data, f, indent=2, sort_keys=True)
