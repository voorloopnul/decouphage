import datetime
import os
import re
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

from src import tools
from src.annotate import Annotate
from src.output import write_gbk


class Pipeline(object):
    def __init__(self, database, input_file, prodigal, threads, output_file):
        self.database = database
        self.contig_file = input_file

        # Operational options
        self.use_prodigal = prodigal
        self.blast_threads = threads
        self.output_file = output_file

        tmp_dir = tempfile.TemporaryDirectory()
        self.tmp_query_file = os.path.join(tmp_dir.name, 'query.fa')
        self.tmp_blast_file = os.path.join(tmp_dir.name, 'blast.tsv')

        self.genome = None
        self.orf_map = None
        self.run()

    def orf_calling(self):
        if self.use_prodigal:
            self.orf_map = tools.run_prodigal(self.contig_file)
        else:
            self.orf_map = tools.run_phanotate(self.contig_file)

    def run(self):
        self.load_genome_from_file()
        self.orf_calling()
        self.load_features()
        self.prepare_query_file()
        tools.run_blast(self.blast_threads, self.tmp_query_file, self.tmp_blast_file, self.database)

        qualifiers = {record["qseqid"]: record for record in Annotate(self.tmp_blast_file).run()}
        self.enrich_features(qualifiers)

        write_gbk(self.genome.values(), self.output_file)

    def enrich_features(self, qualifiers):
        for contig in self.genome.values():
            for feature in contig.features:
                blast_result = qualifiers.get(int(feature.id), {})

                product = blast_result.get("stitle", "Unknown")
                pattern = r"\[.*?\]"
                product = re.sub(pattern, "", product).rstrip()
                product = product.replace("TPA: MAG TPA", "")

                tag = "PREF_{l:04d}".format(l=int(feature.id))

                quals = {
                    "gene": feature.id,
                    "product": product,
                    "locus_tag": tag,
                    "translation": self.cleanup_sequence(feature, contig),
                    "protein_id": blast_result.get("sseqid", "N/A"),
                }
                feature.qualifiers = quals

    def load_features(self):
        for contig_label, orf_list in self.orf_map.items():
            contig = self.genome[contig_label]
            for orf in orf_list:
                idx, start, end, strand = orf.split("_")
                feature = SeqFeature(
                    FeatureLocation(int(start) - 1, int(end)),
                    strand=1 if strand == "+" else -1,
                    type="CDS",
                    qualifiers={},
                    id=idx[1:],
                )
                contig.features.append(feature)

    def load_genome_from_file(self):
        today_date = str(datetime.date.today().strftime("%d-%b-%Y")).upper()
        self.genome = SeqIO.to_dict(SeqIO.parse(self.contig_file, "fasta"))

        for contig in self.genome.values():
            contig.annotations = {"molecule_type": "DNA", "date": today_date}
            contig.id = Path(self.contig_file).stem

    def prepare_query_file(self):
        with open(self.tmp_query_file, "a") as fh:
            for contig in self.genome.values():
                for feature in contig.features:
                    fh.write(f">{feature.id}\n")
                    fh.write(self.cleanup_sequence(feature, contig) + "\n")

    @staticmethod
    def cleanup_sequence(feature, contig):
        sequence = str(feature.extract(contig.seq).translate(table=11))
        if sequence.endswith("*"):
            sequence = sequence[:-1]
        return sequence
