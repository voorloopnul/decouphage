import re

import numpy as np
import pandas as pd
import logging

logger = logging.getLogger(__name__)

header = [
    "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "slen", "stitle", "qlen"
]


class Annotate(object):
    def __init__(self, tmp_blast_file):
        df = pd.read_csv(tmp_blast_file, delimiter="\t", header=None)
        df.columns = header

        # df = df.head(10)

        self.genes_list = list(df["qseqid"].unique())
        self.df = df

    def get_df_for_spcific_gene(self, gene):
        return self.df[self.df["qseqid"] == gene]

    @staticmethod
    def get_df_for_genes_same_length_range(df):
        """
        qlen is the size of the protein we want to annotate.
        slen is the size of the protein in the database.
        """
        offset = int(df["qlen"].unique()[0] * 0.10)

        _df = df[df["slen"] > df["qlen"] - offset]
        _df = _df[_df["slen"] < _df["qlen"] + offset]
        return _df

    @staticmethod
    def remove_species(stitle):
        pattern = r"\[.*?\]"
        return re.sub(pattern, "", stitle).rstrip(" ")

    def default_pass(self, _df, step):
        for index, row in _df.iterrows():
            if row["stitle"] is np.NAN:
                # in case some database entry don't contain a description/product, we skip it.
                logger.warning(f"Missing subject title for blast database entry: {row['sseqid']}")
                continue

            blast_result = row.to_dict()
            blast_result["stitle"] = self.remove_species(blast_result["stitle"])
            logging.debug(f"{step}: {blast_result}")
            return blast_result
        return None

    def run(self):
        qualifiers = []
        for gene in self.genes_list:
            _df_all = self.get_df_for_spcific_gene(gene)
            _df_same_length = self.get_df_for_genes_same_length_range(_df_all)

            blast_result = self.default_pass(_df_same_length, "1st")
            if not blast_result:
                blast_result = self.default_pass(_df_all, "2nd")

            if blast_result:
                qualifiers.append(blast_result)
            else:
                logger.info(f"No good candidate found for CDS #{gene}")

        return qualifiers
