![Tux, the Linux mascot](https://raw.githubusercontent.com/voorloopnul/voorloopnul/357a7ead62584e352c61b008790fe38d4aff5664/logos/decouphage.png)


# Decouphage: the art of decorating a Phage genome by gluing ~~colorful~~ meaningful ~~paper~~ features into it.

## Description

As the name suggests decouphage is a tool designed to annotate phage genomes. It only external dependency is ncbi-blast+
everything else is optional. 
 
### Highlights

 - Can be easily installed in Linux or Mac computers. Only requirement is ncbi-blast+.
 - Offer a second-pass option on top of existing genbank files.
 - Can be extended with prodigal, but as default it uses phanotate for ORF calling. 
 - Decouphage is fast, using a Macbook most phage genomes can be annotated in less than a minute.
 - Uses ncbi NR database containing non-identical sequences from GenBank CDS translations, PDB, Swiss-Prot, PIR, and PRF. 
 
## Why and how should I use decouphage?

1. ### I want to discover and annotate a lot of ORFs:
    
    decouphage genome.fasta -o genome.gb

2. ### I want to use prodigal to find my genes:

    decouphage genome.fasta -o genome.gb --prodigal

3. ### I don't want to call any CDS, I already have it in a genbank and whant to use those:

    > In this mode decouphage will not call any orf, will just generate annotation 

    decouphage genome.gbk -o genome.gb --no-orf-calling

5. ### I have a genbak annotated by other tool and want to fill the blanks!
    
    > In this mode decouphage will not call any orf, it will also preserve the original product and annotate only the
   > ones `Hypothetical Protein` and `Phage protein`


    decouphage genome.gbk -o genome.gb --second-pass


Doens't matter what annotation tool you choose to use, if you are manually investigating a genome you are going to try 
and blast all those `Hypotethical Protein` or `Phage protein` that are missing in your annotation.




3. ###

## Installation

You have multiple options to install and run decouphage:

### Ubuntu

Install decouphage:

    pip install decouphage

Install ncbi-blast+
    
    apt install ncbi-blast+

Optional: Install dependencies:

    apt install prodigal

### Docker

Run with docker (Already includes dependencies and databases):

    docker run decouphage/decouphage

## Databases

Decouphage database is derived from NCBI NR database clustered at 90% identity and 90% sequence length.

### Downloading database

Download database to default location in $HOME/.decouphage/db/

    decouphage --download_db

Download database to a different path:

    decouphage --download_db /some/location/to/store/database 

### Making custom databases

Make blast database

    makeblastdb -in database.fa -parse_seqids -blastdb_version 5 -dbtype prot

# Usage

## Options

    Usage: decouphage [OPTIONS] INPUT_FILE
    
    Options:
      --prodigal             Use prodigal for orf calling instead of phanotate.
      -d, --db PATH
      -o, --output TEXT
      -t, --threads INTEGER  [default: 1]
      --tmp_dir TEXT         Folder for intermediate files.
      --no-orf-calling       Annotate CDS from genbank file.
      --help                 Show this message and exit.


## Examples

Using prodigal:

    $ decouphage --db nr --prodigal --threads 4 input_genome.fa

Using phanotate (default):

    $ decouphage --db nr --threads 4 genome.fasta
