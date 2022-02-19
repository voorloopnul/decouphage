import re

import pandas as pd

from src.output import write_demo

header = [
    "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "slen", "stitle", "qlen"
]


class Annotate(object):
    def __init__(self):
        df = pd.read_csv("tmp/blast.tsv", delimiter="\t", header=None)
        df.columns = header

        # df = df.head(10)
        df[["gene", "start", "end", "strand"]] = df["qseqid"].str.split(pat="_", expand=True)
        _start = df["start"].astype(int)
        _end = df["end"].astype(int)
        query_len = _end - (_start - 1)  # adjust for python 0-index
        query_len = query_len - 3  # subtract stop codon
        query_len = query_len / 3  # get protein query size
        df["query_len"] = query_len.astype(int)

        #del df["start"]
        #del df["end"]

        self.genes_list = list(df["gene"].unique())
        self.df = df

    def get_df_for_spcific_gene(self, gene):
        return self.df[self.df["gene"] == gene]

    @staticmethod
    def get_df_for_genes_same_length_range(df):
        """
        query_len is the size of the protein we want to annotate.
        slen is the size of the protein in the database.
        """
        offset = int(df["query_len"].unique()[0] * 0.10)

        _df = df[df["slen"] > df["query_len"] - offset]
        _df = _df[_df["slen"] < _df["query_len"] + offset]
        return _df

    def run(self):
        qualifiers = []
        for gene in self.genes_list:
            _df = self.get_df_for_spcific_gene(gene)
            _df = self.get_df_for_genes_same_length_range(_df)
            #print(_df)

            if _df.empty:
                print(f"{gene} Empty...")
                #_df = self.get_df_for_spcific_gene(gene)

            for index, row in _df.iterrows():
                if "hypothetical" not in row['stitle']:
                    if "putative" not in row["stitle"]:
                        print(row["gene"], row["query_len"], row["slen"], row["stitle"])
                        blast_result = row.to_dict()
                        blast_result["stitle"] = re.sub("[\(\[].*?[\)\]]", "", blast_result["stitle"])
                        qualifiers.append(blast_result)
                        break
            else:
                for index, row in _df.iterrows():
                    print(row["gene"], row["query_len"], row["slen"], row["stitle"], "LOW")
                    blast_result = row.to_dict()
                    blast_result["stitle"] = re.sub("[\(\[].*?[\)\]]", "", blast_result["stitle"])
                    qualifiers.append(blast_result)
                    break
        return qualifiers

if __name__ == "__main__":
    Annotate().run()
