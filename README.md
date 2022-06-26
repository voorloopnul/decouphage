![Tux, the Linux mascot](https://raw.githubusercontent.com/voorloopnul/voorloopnul/357a7ead62584e352c61b008790fe38d4aff5664/logos/decouphage.png)


# Decouphage: the art of decorating a Phage genome by gluing ~~colorful~~ meaningful feature cutouts into it.

## Description

As the name suggests decouphage is a tool designed to annotate phage genomes. It only external dependency is ncbi-blast+
everything else is optional. 
 
### Highlights

 - Can be easily installed in Linux or Mac computers. Only requirement is ncbi-blast+.
 - Can be extended with prodigal
 - Decouphage is fast, using a Macbook most phage genomes can be annotated in less than 5 minutes.
 - Uses ncbi NR database containing non-identical sequences from GenBank CDS translations, PDB, Swiss-Prot, PIR, and PRF. 
 
## Installation

You have multiple options to install and run decouphage:

### Ubuntu

Install decouphage:

    pip install decouphage

Install ncbi-blast+
    
    apt install ncbi-blast+ prodigal

Optional: Install dependencies:

    apt install prodigal

### Docker

Run with docker (Already includes dependencies and databases):

    docker run decouphage/decouphage

## Databases

Decouphage database is derived from NCBI NR.

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
