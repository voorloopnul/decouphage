

def write_gbk(genome):
    for contig in genome:
        print(contig.format("gb"))
        pass

    #SeqIO.write(rec, handle, "genbank")
    # handle.close()

