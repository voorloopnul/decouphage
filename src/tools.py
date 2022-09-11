import subprocess
from collections import defaultdict
import logging
import os
import sys

BLAST_FMT = "qseqid sseqid pident length evalue bitscore slen stitle qlen"
logger = logging.getLogger(__name__)


def run_blast(threads, query_file, blast_file, database):
    if not os.path.exists(database):
        logging.error("Default database not found.")
        sys.exit(1)

    blast_cmd = [
        "blastp",
        "-db", database,
        "-query", query_file,
        "-evalue", "1e-5",
        "-outfmt", f"'6 {BLAST_FMT}'",
        "-num_threads", f"{threads}",
        "-task", "blastp-fast",
        "-out", blast_file
    ]

    try:
        rt = subprocess.Popen(" ".join(blast_cmd), shell=True)
        rt.communicate()

    except FileNotFoundError:
        logging.error("blast not found!")
        sys.exit(1)


def run_prodigal(path):
    """
    Expects as input a FASTA file and return the prodigal default output as a list:
    ['>1_70_930_+', '>2_1282_2259_+', '>3_2276_2689_+', '>4_3116_3667_-', ... ]
    """
    prodigal_cmd = f"prodigal -i {path} -f sco -p meta"
    try:
        rt = subprocess.Popen(prodigal_cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        logging.error("Prodigal not found!")
        sys.exit(1)

    output, error = rt.communicate()
    output = output.decode("utf-8")
    features = defaultdict(list)

    idx = 1
    for line in output.split("\n"):
        if "seqhdr" in line:
            contig = line.split('"')[1].split(" ")[0]
        if line.startswith(">"):
            orf = line.split("_")
            orf[0] = f">{idx}"
            idx += 1
            features[contig].append("_".join(orf))

    return features


def run_phanotate(path):
    phanotate_cmd = f"phanotate.py {path}"
    rt = os.popen(phanotate_cmd).read()
    rt = rt.split("\n")

    features = defaultdict(list)
    idx = 1
    for line in rt:
        if line.split("\t") == [""]:
            break

        if not line.startswith("#"):
            orf = line.split("\t")[0:4]
            strand = orf[2]
            contig = orf[3]
            start, end = sorted([int(x) for x in orf[0:2]])
            orf = "_".join([f">{idx}", str(start), str(end), strand])
            idx += 1
            features[contig].append(orf)

    return features


def run_trnascan(path):
    """
    """

    trnascan_cmd = f"tRNAscan-SE -B -q --brief {path}"
    rt = os.popen(trnascan_cmd).read()
    rt = rt.split("\n")

    features = defaultdict(list)

    for line in rt[:-1]:
        trna = line.split("\t")
        orf = "_".join([f">{trna[1]}", str(trna[2].strip()), str(trna[3].strip()), trna[4].strip(), trna[5].strip()])
        features[trna[0].strip()].append(orf)

    return features
