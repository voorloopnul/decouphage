import pandas as pd

header = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue",
    "bitscore", "slen", "stitle"
]


class Annotate(object):
    def __init__(self):
        df = pd.read_csv("tmp/blast.tsv", delimiter="\t", header=None)
        df.columns = header

        # df = df.head(10)
        df[["gene", "start", "end", "strand"]] = df["qseqid"].str.split(pat="_", expand=True)
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)
        df["query_len"] = df["end"] - df["start"]
        df["query_len"] = (df["query_len"] / 3).astype(int)

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
        for gene in self.genes_list:
            _df = self.get_df_for_spcific_gene(gene)
            _df = self.get_df_for_genes_same_length_range(_df)

            if _df.empty:
                print("Empty...")
                _df = self.get_df_for_spcific_gene(gene)

            for index, row in _df.iterrows():
                if "hypothetical" not in row['stitle']:
                    if "putative" not in row["stitle"]:
                        # print(_df.loc[[index]])
                        print(row["gene"], row["query_len"], row["slen"], row["stitle"])
                        break
            else:
                for index, row in _df.iterrows():
                    print(row["gene"], row["query_len"], row["slen"], row["stitle"], "LOW")
                    break


if __name__ == "__main__":
    Annotate().run()
