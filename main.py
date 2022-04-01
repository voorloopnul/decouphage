import os
import sys
from Bio import SeqIO
from src.annotate import Annotate
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re
from src import tools
from src.output import write_gbk


class Pipeline(object):
    def __init__(self, contig_file):
        self.contig_file = contig_file
        self.genome = None
        self.orf_map = None
        self.run()

    def run(self):
        self.load_genome()
        self.cleanup_query_file()
        self.orf_map = tools.run_prodigal(self.contig_file)
        self.load_features()
        self.prepare_query_file()
        tools.run_blast()

        qualifiers = Annotate().run()
        qualifiers = {record["qseqid"]: record for record in qualifiers}

        self.enrich_features(qualifiers)
        write_gbk(self.genome.values())

    def enrich_features(self, qualifiers):
        for contig in self.genome.values():
            for feature in contig.features:
                blast_result = qualifiers.get(int(feature.id), {})

                product = blast_result.get("stitle", "Unknown")
                pattern = r'\[.*?\]'
                product = re.sub(pattern, '', product).rstrip()
                product = product.replace("TPA: MAG TPA", "")

                tag = 'PREF_{l:04d}'.format(l=int(feature.id))

                quals = {
                    'gene': feature.id,
                    'product': product,
                    'locus_tag': tag,
                    'translation': str(feature.extract(contig.seq).translate()),
                    "protein_id": blast_result.get("sseqid", "N/A")
                }
                feature.qualifiers = quals

    def load_features(self):
        for contig_label, orf_list in self.orf_map.items():
            print(contig_label)
            contig = self.genome[contig_label]
            for orf in orf_list:
                idx, start, end, strand = orf.split("_")
                feature = SeqFeature(
                    FeatureLocation(int(start) - 1, int(end)),
                    strand=1 if strand == "+" else -1,
                    type="CDS",
                    qualifiers={},
                    id=idx[1:]
                )
                contig.features.append(feature)

    def load_genome(self):
        self.genome = SeqIO.to_dict(SeqIO.parse(self.contig_file, "fasta"))
        for contig in self.genome.values():
            contig.annotations = {"molecule_type": "DNA"}

    def prepare_query_file(self):
        with open("tmp/query.fa", "a") as fh:
            for contig in self.genome.values():
                for feature in contig.features:
                    fh.write(f">{feature.id}\n")
                    fh.write(str(feature.extract(contig.seq).translate()) + "\n")

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
    input_path = sys.argv[1]
    Pipeline(input_path)
