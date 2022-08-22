from Bio import SeqIO
import plotly.graph_objects as go

word_of_interest = [
    "endonuclease",
    "exonuclease",
    "helicase",
    "hydrolase",
    "kinase",
    "ligase",
    "methyltransferase",
    "polymerase",
    "primase",
    "protease",
    "recombinase",
    "reductase",
    "synthase",
    "terminase",
    "transferase",
    "hypothetical protein",
]


def get_list_of_products(genbank_file):
    product_list = []
    with open(genbank_file) as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    product_list.append(feature.qualifiers.get("product")[0])
    print(f"Found {len(product_list)} products in {genbank_file}")
    return product_list


def count_enzyme_occurrences(decouphage_list, rast_list):
    d_list, r_list, a_list = [], [], []
    for word in word_of_interest:
        d_count, r_count, a_count = 0, 0, 0

        for d, r in zip(decouphage_list, rast_list):
            # count enzyme occurrence in Decouphage
            d_count += 1 if word in d.lower() else 0
            # count enzyme occurrence in rast
            r_count += 1 if word in r.lower() else 0
            # count enzyme occurrence agreement
            a_count += 1 if word in d.lower() and word in r.lower() else 0

        d_list.append(d_count)
        r_list.append(r_count)
        a_list.append(a_count)
    return d_list, r_list, a_list


def generate_markdown_table():
    print("| Enzyme | Decouphage | Rast | Decouphage agreement rate |")
    print("| ------ | ---------- | ---- | ------------------------- |")
    for w, d, r, a in zip(word_of_interest, decouphage_counts, rast_counts, agreement_counts):
        print(f"| {w} | {d:4d} | {r:4d} | {int(a/r*100)}% |")


def generate_figure_01():
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=word_of_interest[:-1],
        y=decouphage_counts[:-1],
        text=decouphage_counts[:-1],
        name='Decouphage',
    ))
    fig.add_trace(go.Bar(
        x=word_of_interest[:-1],
        y=rast_counts[:-1],
        text=rast_counts[:-1],
        name='RAST',
    ))

    # Here we modify the tickangle of the xaxis, resulting in rotated labels.
    fig.update_layout(
        barmode='group',
        xaxis_tickangle=-45,
        font=dict(
            family="Courier New, monospace",
            size=18,
        ),
        autosize=False,
        width=1200,
        height=600,
        title="Product annotation counts"
    )

    fig.show()


product_list_rast = get_list_of_products("100_genomes_rast.gbk")
product_list_decouphage = get_list_of_products("100_genomes_decouphage.gbk")
decouphage_counts, rast_counts, agreement_counts = count_enzyme_occurrences(product_list_decouphage, product_list_rast)

generate_markdown_table()
generate_figure_01()
