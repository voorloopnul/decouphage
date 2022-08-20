# Validation steps

1. **Download genomes from inphared**

    
    wget http://inphared.s3.climb.ac.uk/1Aug2022_genomes_excluding_refseq.fa

2. **Take a sample of 100 genomes**

    
    seqtk sample 1Aug2022_genomes_excluding_refseq.fa 100 > 100_genomes.fa

3. **Annotate it using RAST web interface**


    https://rast.nmpdr.org/rast.cgi

4. **Run the output with decouphage using the --no_orf_calling flag**

    
    decouphage 100_genomes_rast.gbk --no-orf-calling -t 8 -o 100_genomes_decouphage.gbk

5. **Compare the results**

