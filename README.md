![Decouphage logo](https://raw.githubusercontent.com/voorloopnul/voorloopnul/357a7ead62584e352c61b008790fe38d4aff5664/logos/decouphage.png)


# Decouphage: the art of decorating a Phage genome by gluing feature cutouts into it.

### Relevant branches

 - [main branch](https://github.com/voorloopnul/decouphage/tree/main): stable version available in [pypi](https://pypi.org/project/decouphage/) and [dockerhub](https://hub.docker.com/r/voorloop/decouphage).
 - [dev branch](https://github.com/voorloopnul/decouphage/tree/dev): development branch with most recent version of decouphage.
 

## Description

As the name suggests decouphage is a tool designed to annotate phage genomes. It only external dependency is ncbi-blast+
everything else is optional. 
 
### Highlights

 - Can be easily installed in Linux or Mac computers. Only requirement is ncbi-blast+.
 - Can be extended with prodigal, but as default it uses phanotate for ORF calling. 
 - Decouphage is fast, using a Macbook **most phage genomes can be annotated in less than a minute**.
 - Uses ncbi NR database containing non-identical sequences from GenBank CDS translations, PDB, Swiss-Prot, PIR, and PRF. 
 
## How can I use decouphage?

### Options

    Usage: decouphage [OPTIONS] INPUT_FILE
    
    Options:
      --prodigal             Use prodigal for orf calling instead of phanotate.
      -d, --db PATH
      -o, --output TEXT
      -t, --threads INTEGER  [default: 1]
      --tmp_dir TEXT         Folder for intermediate files.
      --no-orf-calling       Annotate CDS from genbank file.
      --help                 Show this message and exit.


### I want to discover and annotate a lot of ORFs:
 
    decouphage genome.fasta -o genome.gb

### I want to use prodigal to find my genes:

    decouphage genome.fasta -o genome.gb --prodigal

### I have a genbank with poor annotation and want more:

In this mode decouphage will reuse the genbank ORFs and just run the annotation procedure.

    decouphage genome.gbk -o genome.gb --no-orf-calling

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
