from Bio import SeqIO

word_of_interest = [
    "baseplate",
    "endonuclease",
    "exonuclease",
    "head",
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
    "thioredoxin",
    "tail",
#    "hypothetical protein",
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


product_list_rast = get_list_of_products("100_genomes_rast.gbk")
product_list_decouphage = get_list_of_products("100_genomes_decouphage.gbk")

r_s, d_s = 0, 0
for word in word_of_interest:
    r_count, d_count, a_count = 0, 0, 0
    for r, d in zip(product_list_rast, product_list_decouphage):
        r_count += 1 if word in r.lower() else 0
        d_count += 1 if word in d.lower() else 0
        a_count += 1 if word in d.lower() and word in r.lower() else 0
    print(f"R: {r_count:4d} | D: {d_count:4d} | DaR: {int(a_count/r_count*100)}% | {word}")
    r_s += r_count
    d_s += d_count


print(f"R: {r_s:4d} | D: {d_s:4d} ")


