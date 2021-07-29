import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from src.annotate import Annotate

BLAST_FMT = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen stitle"


def get_description(uuid):
    desc_cmd = f"blastdbcmd -db db/database.fa -entry '{uuid}'"
    entry = os.popen(desc_cmd).read()
    header = entry.split("\n")[0]
    return ' '.join(header.split(" ")[1:3])


def run_blast():
    print("blasting...")
    blast_cmd = f"blastp -db db/database.fa -query tmp/query.fa -evalue 1e-5 -outfmt '6 {BLAST_FMT}' -num_threads 8 -out tmp/blast.tsv"
    os.system(blast_cmd)


def run_prodigal(path):
    prodigal_cmd = f"prodigal -i {path} -p meta -f sco"
    rt = os.popen(prodigal_cmd).read()
    rt = rt.split("\n")
    return rt[2:-1]


class Pipeline(object):
    def __init__(self, contig_file):
        self.contig_file = contig_file
        self.genome = None
        self.genes = None

    def run(self):
        self.load_genome()
        self.cleanup_query_file()
        self.genes = run_prodigal(self.contig_file)
        self.prepare_query_file()
        run_blast()
        Annotate().run()

    def load_genome(self):
        print(self.contig_file)
        for record in SeqIO.parse(self.contig_file, "fasta"):
            self.genome = record
            break

    def prepare_query_file(self):
        for gene in self.genes:
            # print(gene)
            pos, start, end, strand = gene.split("_")
            fixed_start = int(start) - 1  # python is 0 index based
            fixed_end = int(end)

            sequence = self.genome.seq[fixed_start:fixed_end]
            sequence = Seq(sequence)
            if strand == "-":
                sequence = sequence.reverse_complement()

            # print(sequence.translate())

            with open("tmp/query.fa", "a") as fh:
                fh.write(f"{gene}\n")
                fh.write(str(sequence.translate()[0:-1]) + "\n")

    @staticmethod
    def cleanup_query_file():
        try:
            os.mkdir("tmp")
        except Exception as err:
            pass

        with open("tmp/query.fa", "w") as fh:
            fh.write("")


if __name__ == '__main__':
    print("Decouphage 0.1")
    contigs_path = sys.argv[1]

    Pipeline(contigs_path).run()
