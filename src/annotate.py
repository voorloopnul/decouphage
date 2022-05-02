import re
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

header = [
    "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "slen", "stitle", "qlen"
]


class Annotate(object):
    def __init__(self):
        df = pd.read_csv("tmp/blast.tsv", delimiter="\t", header=None)
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
        return re.sub("[\(\[].*?[\)\]]", "", stitle)

    def first_pass(self, _df):
        bad_words = ["hypothetical", "putative", "unknown"]
        for index, row in _df.iterrows():
            if not any(word in row["stitle"] for word in bad_words):
                blast_result = row.to_dict()
                blast_result["stitle"] = self.remove_species(blast_result["stitle"])
                logging.debug(f"1st: {blast_result}")
                return blast_result
        return None

    def second_pass(self, _df):
        for index, row in _df.iterrows():
            blast_result = row.to_dict()
            blast_result["stitle"] = self.remove_species(blast_result["stitle"])
            logging.debug(f"2nd: {blast_result}")
            return blast_result
        return None

    def third_pass(self, _df):
        bad_words = ["hypothetical", "putative", "unknown"]
        for index, row in _df.iterrows():
            if not any(word in row["stitle"] for word in bad_words):
                blast_result = row.to_dict()
                blast_result["stitle"] = self.remove_species(blast_result["stitle"])
                logging.debug(f"3rd: {blast_result}")
                return blast_result
        return None

    def run(self):
        qualifiers = []
        for gene in self.genes_list:
            _df = self.get_df_for_spcific_gene(gene)
            _df = self.get_df_for_genes_same_length_range(_df)

            blast_result = self.first_pass(_df)

            if not blast_result:
                blast_result = self.second_pass(_df)

            if not blast_result:
                blast_result = self.third_pass(self.get_df_for_spcific_gene(gene))

            if blast_result:
                qualifiers.append(blast_result)
            else:
                logger.info(f"No good candidate found for CDS #{gene}")


        return qualifiers


if __name__ == "__main__":
    Annotate().run()
