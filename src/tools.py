from collections import defaultdict
import logging
import os

logging.basicConfig(level=logging.INFO)
BLAST_FMT = "qseqid sseqid pident length evalue bitscore slen stitle qlen"
logger = logging.getLogger(__name__)


def run_blast(threads, query_file, blast_file):
    logger.info("Blasting sequences...")

    blast_cmd = [
        "blastp",
        "-db", "db/nr.fa",
        "-query", query_file,
        "-evalue", "1e-5",
        "-outfmt", f"'6 {BLAST_FMT}'",
        "-num_threads", f"{threads}",
        "-out", blast_file
    ]

    os.system(" ".join(blast_cmd))


def run_prodigal(path):
    """
    Expects as input a FASTA file and return the prodigal default output as a list:
    ['>1_70_930_+', '>2_1282_2259_+', '>3_2276_2689_+', '>4_3116_3667_-', ... ]
    """
    prodigal_cmd = f"prodigal -i {path} -f sco -p meta"
    rt = os.popen(prodigal_cmd).read()
    rt = rt.split("\n")

    features = defaultdict(list)
    idx = 1
    for line in rt:
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
