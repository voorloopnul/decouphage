import re
import pandas as pd

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

    def run(self):
        qualifiers = []
        for gene in self.genes_list:
            _df = self.get_df_for_spcific_gene(gene)
            _df = self.get_df_for_genes_same_length_range(_df)

            if _df.empty:
                print(f"{gene} Not a good candidate found...")
                #_df = self.get_df_for_spcific_gene(gene)

            for index, row in _df.iterrows():
                if "hypothetical" not in row['stitle']:
                    if "putative" not in row["stitle"]:
                        if "unknown" not in row["stitle"]:
                            print("A", row.to_dict())
                            blast_result = row.to_dict()
                            blast_result["stitle"] = re.sub("[\(\[].*?[\)\]]", "", blast_result["stitle"])
                            qualifiers.append(blast_result)
                            break
            else:
                for index, row in _df.iterrows():
                    print("B", row.to_dict())
                    blast_result = row.to_dict()
                    blast_result["stitle"] = re.sub("[\(\[].*?[\)\]]", "", blast_result["stitle"])
                    qualifiers.append(blast_result)
                    break
        return qualifiers


if __name__ == "__main__":
    Annotate().run()
