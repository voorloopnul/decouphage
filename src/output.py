from Bio import SeqIO


def write_gbk(genome, output_file):
    with open(output_file, "w") as fh:
        for contig in genome:
            SeqIO.write(contig, fh, "genbank")

