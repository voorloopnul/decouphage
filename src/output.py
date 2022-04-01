from Bio import SeqIO


def write_gbk(genome):
    with open("output.gbk", "w") as fh:
        for contig in genome:
            SeqIO.write(contig, fh, "genbank")

