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


    ./python compare.py


## Example output

![alt text](../assets/decouphage_image_01.png?raw=true)


| Enzyme | Decouphage | Rast | Decouphage agreement rate |
| ------ | ---------- | ---- | ------------------------- |
| endonuclease |  232 |   68 | 94% |
| exonuclease |   59 |   29 | 58% |
| helicase |  113 |   71 | 70% |
| hydrolase |   90 |   26 | 73% |
| kinase |   66 |   23 | 86% |
| ligase |   52 |   19 | 94% |
| methyltransferase |   54 |   23 | 65% |
| polymerase |  153 |   80 | 76% |
| primase |   64 |   37 | 78% |
| protease |   68 |   20 | 85% |
| recombinase |   36 |    7 | 28% |
| reductase |   94 |   65 | 90% |
| synthase |   43 |   33 | 84% |
| terminase |  139 |   92 | 94% |
| transferase |  159 |   43 | 60% |
| hypothetical protein | 3945 | 6302 | 53% |
| phage protein |    0 | 1626 | 0% |
