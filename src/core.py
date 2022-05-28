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


def probe_filetype(input_file):
    with open(input_file, 'rb') as fh:
        data = fh.read(2)

    if data.startswith(b">"):
        return "FASTA"
    else:
        return "GENBANK"


def get_database_default_path():
    return Path.home() / Path(".decouphage/db/nr/nr.fa")
