import datetime
import logging
import os
import re
import shutil
import sys
import tempfile
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from src import tools
from src.annotate import Annotate
from src.core import cds_calling_from_genbank, get_database_default_path, validate_input
from src.output import write_gbk

logger = logging.getLogger(__name__)
total_steps = 9


class Pipeline(object):
    def __init__(self, database: str, input_file: str, prodigal: bool, threads: int, output_file: str, tmp_dir: str,
                 merge_gbk: bool, locus_tag: str):

        self.database = database if database else str(get_database_default_path())
        self.contig_file = input_file
        self.locus_tag = locus_tag

        # Operational options
        self.use_prodigal = prodigal
        self.blast_threads = threads
        self.output_file = output_file
        self.merge_gbk = merge_gbk

        tmp_dir_handler = tempfile.TemporaryDirectory()
        self.tmp_dir = tmp_dir if tmp_dir else tmp_dir_handler.name
        if tmp_dir:
            os.mkdir(tmp_dir)

        self.tmp_input_gbk_file = os.path.join(self.tmp_dir, 'input.gbk')  # genbank
        self.tmp_input_fna_file = os.path.join(self.tmp_dir, 'input.fna')  # contigs
        self.tmp_query_faa_file = os.path.join(self.tmp_dir, 'query.faa')  # ORFs
        self.tmp_blast_tsv_file = os.path.join(self.tmp_dir, 'blast.tsv')  # Blast result

        self.genome = None
        self.orf_map = None
        self.trna_map = None
        self.run()

    def run(self):
        self.load_genome_from_file()
        self.add_contig_metadata()

        logger.info(f"[2/{total_steps}] Starting ORF calling")
        self.orf_map = self.orf_calling()

        logger.info(f"[3/{total_steps}] Starting tRNA calling")
        self.trna_map = self.rna_calling()

        logger.info(f"[4/{total_steps}] Loading features into biopython SeqFeature")
        self.load_features()

        logger.info(f"[5/{total_steps}] Preparing blast query file")
        self.prepare_query_file()

        logger.info(f"[6/{total_steps}] Running blast with {self.database}")
        tools.run_blast(self.blast_threads, self.tmp_query_faa_file, self.tmp_blast_tsv_file, self.database)

        logger.info(f"[7/{total_steps}] Selecting best hits")
        qualifiers = {record["qseqid"]: record for record in Annotate(self.tmp_blast_tsv_file).run()}

        logger.info(f"[8/{total_steps}] Enriching features")
        self.enrich_features(qualifiers)

        logger.info(f"[9/{total_steps}] Writing output.")
        write_gbk(self.genome.values(), self.output_file)

    def load_genome_from_file(self):
        input_type = validate_input(self.contig_file)
        if not input_type:
            logger.error("Failed to detect input type")
            sys.exit(1)

        logger.info(f"[1/{total_steps}] Input type inferred as {input_type}")

        if input_type == 'FASTA' and self.merge_gbk:
            logger.error("Can't disable ORF calling if input is a FASTA file.")
            sys.exit(1)

        if input_type == 'FASTA':
            shutil.copyfile(self.contig_file, self.tmp_input_fna_file)
        elif input_type == 'GENBANK':
            with open(self.contig_file, "r") as input_handle:
                sequences = SeqIO.parse(input_handle, "genbank")
                with open(self.tmp_input_fna_file, "w") as output_handle:
                    SeqIO.write(sequences, output_handle, "fasta")
            shutil.copyfile(self.contig_file, self.tmp_input_gbk_file)
        else:
            logger.error("The code should never have hit this else!")
            sys.exit(1)

        self.contig_file = self.tmp_input_fna_file
        try:
            self.genome = SeqIO.to_dict(SeqIO.parse(self.contig_file, "fasta"))
        except ValueError as e:
            duplicated_contig = re.findall(r"'(.*?)'", str(e), re.DOTALL)[0]
            logger.error(f"A duplicated record was found: {duplicated_contig}")
            sys.exit(1)

    def add_contig_metadata(self):
        today_date = str(datetime.date.today().strftime("%d-%b-%Y")).upper()
        for contig in self.genome.values():
            contig.annotations = {"molecule_type": "DNA", "date": today_date, "accessions": ""}
            # Leave .id empty to mimic prokka output.
            contig.id = ""

    def rna_calling(self):
        return tools.run_trnascan(self.tmp_input_fna_file)

    def orf_calling(self):
        # Don't call ORFs, just copy it from the input gbk file.
        if self.merge_gbk:
            logger.warning(f" - decouphage will skip this step and use CDS from Genbank file.")
            return cds_calling_from_genbank(self.tmp_input_gbk_file)

        # Use prodigal for ORF calling
        if self.use_prodigal:
            logger.info(f" - with Prodigal")
            return tools.run_prodigal(self.tmp_input_fna_file)

        # Use phanotate for ORF calling
        logger.info(f" - with Phanotate")
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
                )
                contig.features.append(feature)

            for trna in self.trna_map[contig_label]:
                idx, start, end, rna_type, anti_codon = trna.split("_")
                start, end = int(start), int(end)
                strand = 1 if end > start else -1
                positions = [start, end]
                feature = SeqFeature(
                    FeatureLocation(min(positions) - 1, max(positions)),
                    strand=strand,
                    type="tRNA",
                    qualifiers={"product": f"{rna_type}-{anti_codon}"},
                )
                contig.features.append(feature)

        # order features based on start position
        for contig in self.genome.values():
            contig.features.sort(key=lambda x: x.location.start, reverse=False)

        # give feature id an incremental count
        idx = 1
        for contig in self.genome.values():
            for feature in contig.features:
                feature.id = idx
                idx += 1

    def prepare_query_file(self):
        with open(self.tmp_query_faa_file, "a") as fh:
            for contig in self.genome.values():
                for feature in contig.features:
                    if feature.type == "CDS":
                        fh.write(f">{feature.id}\n")
                        fh.write(self.clean_sequence(feature, contig) + "\n")

    def enrich_features(self, qualifiers):
        for contig in self.genome.values():
            for feature in contig.features:
                tag = f"{self.locus_tag}_{int(feature.id):04d}"
                quals = feature.qualifiers
                quals["feature"] = feature.id
                quals["locus_tag"] = tag

                if feature.type == "CDS":
                    blast_result = qualifiers.get(int(feature.id), {})
                    quals["product"] = blast_result.get("stitle", "hypothetical protein"),
                    quals["translation"] = self.clean_sequence(feature, contig)
                    quals["protein_id"] = blast_result.get("sseqid", "N/A")
                feature.qualifiers = quals

    @staticmethod
    def clean_sequence(feature, contig):
        sequence = str(feature.extract(contig.seq).translate(table=11))
        if sequence.endswith("*"):
            sequence = sequence[:-1]
        return sequence
