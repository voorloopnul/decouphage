from collections import defaultdict
from pathlib import Path
from Bio import SeqIO


def cds_calling_from_genbank(input_handle):
    features = defaultdict(list)

    contigs = SeqIO.parse(input_handle, "genbank")
    idx = 1
    for contig in contigs:
        for feature in contig.features:
            if feature.type == "CDS":
                strand = "+" if feature.location.strand == 1 else "-"
                orf = f">{idx}_{feature.location.start+1}_{feature.location.end}_{strand}"
                features[contig.id].append(orf)
                idx += 1
    return features


def get_database_default_path():
    return Path.home() / Path(".decouphage/db/ncbi_phages/ncbi_phages.fa")


def is_fasta(file_path):
    has_sequence = False

    fasta = SeqIO.parse(file_path, "fasta")
    for i in fasta:
        if len(i.seq) > 0:
            has_sequence = True
            return has_sequence
    return has_sequence


def is_genbank(file_path):
    has_sequence = False
    fasta = SeqIO.parse(file_path, "genbank")
    for i in fasta:
        if len(i.seq) > 0:
            has_sequence = True
            return has_sequence
    return has_sequence


def validate_input(input_file):
    try:
        if is_fasta(input_file):
            return "FASTA"
    except UnicodeDecodeError:
        return None

    try:
        if is_genbank(input_file):
            return "GENBANK"
    except UnicodeDecodeError:
        return None
