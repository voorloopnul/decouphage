import datetime
import os
import re
import shutil
import tempfile
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from src import tools
from src.annotate import Annotate
from src.core import cds_calling_from_genbank, probe_filetype
from src.output import write_gbk


class Pipeline(object):
    def __init__(self, database: str, input_file: str, prodigal: bool, threads: int, output_file: str, tmp_dir: str,
                 merge_gbk: bool):

        self.database = database if database else "db/nr/nr.fa"
        self.contig_file = input_file

        # Operational options
        self.use_prodigal = prodigal
        self.blast_threads = threads
        self.output_file = output_file
        self.merge_gbk = merge_gbk

        tmp_dir_handler = tempfile.TemporaryDirectory()
        self.tmp_dir = tmp_dir if tmp_dir else tmp_dir_handler.name

        self.tmp_input_gbk_file = os.path.join(self.tmp_dir, 'input.gbk')  # genbank
        self.tmp_input_fna_file = os.path.join(self.tmp_dir, 'input.fna')  # contigs
        self.tmp_query_faa_file = os.path.join(self.tmp_dir, 'query.faa')  # ORFs
        self.tmp_blast_tsv_file = os.path.join(self.tmp_dir, 'blast.tsv')  # Blast result

        self.genome = None
        self.orf_map = None
        self.run()

    def run(self):
        self.load_genome_from_file()
        self.annotate_contigs()
        self.orf_map = self.orf_calling()
        self.load_features()
        self.prepare_query_file()
        tools.run_blast(self.blast_threads, self.tmp_query_faa_file, self.tmp_blast_tsv_file, self.database)

        qualifiers = {record["qseqid"]: record for record in Annotate(self.tmp_blast_tsv_file).run()}
        self.enrich_features(qualifiers)

        write_gbk(self.genome.values(), self.output_file)

    def load_genome_from_file(self):
        input_type = probe_filetype(self.contig_file)

        if input_type == 'FASTA':
            shutil.copyfile(self.contig_file, self.tmp_input_fna_file)
        else:
            with open(self.contig_file, "r") as input_handle:
                sequences = SeqIO.parse(input_handle, "genbank")
                with open(self.tmp_input_fna_file, "w") as output_handle:
                    SeqIO.write(sequences, output_handle, "fasta")
            shutil.copyfile(self.contig_file, self.tmp_input_gbk_file)

        self.contig_file = self.tmp_input_fna_file
        self.genome = SeqIO.to_dict(SeqIO.parse(self.contig_file, "fasta"))

    def annotate_contigs(self):
        today_date = str(datetime.date.today().strftime("%d-%b-%Y")).upper()
        for contig in self.genome.values():
            contig.annotations = {"molecule_type": "DNA", "date": today_date}
            contig.id = Path(self.contig_file).stem

    def orf_calling(self):
        # Don't call ORFs, just copy it from the input gbk file.
        if self.merge_gbk:
            return cds_calling_from_genbank(self.tmp_input_gbk_file)

        # Use prodigal for ORF calling
        if self.use_prodigal:
            return tools.run_prodigal(self.tmp_input_fna_file)

        # Use phanotate for ORF calling
        return tools.run_phanotate(self.tmp_input_fna_file)

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

    def prepare_query_file(self):
        with open(self.tmp_query_faa_file, "a") as fh:
            for contig in self.genome.values():
                for feature in contig.features:
                    fh.write(f">{feature.id}\n")
                    fh.write(self.clean_sequence(feature, contig) + "\n")

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
                    "translation": self.clean_sequence(feature, contig),
                    "protein_id": blast_result.get("sseqid", "N/A"),
                }
                feature.qualifiers = quals

    @staticmethod
    def clean_sequence(feature, contig):
        sequence = str(feature.extract(contig.seq).translate(table=11))
        if sequence.endswith("*"):
            sequence = sequence[:-1]
        return sequence
