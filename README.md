# Decouphage: the art of decorating a Phage genome by gluing colorful features into it.

## Description

As the name suggests decouphage is a tool designed to annotate phage genomes. This decision has been made
in order to tweak the annotation process regarding tools, preferences & databases and to streamline further development &
maintenance of the software.
 
### Highlights

 - Standalone tool that can easily be installed in Linux or Mac computers, can be extended with extra tools like blast
and prodigal.
 - Decouphage is fast, using blast and prodigal it will annotate most phage genomes in less than 5 minutes, using the 
embedded aligner takes longer but still below 10 minutes per genome.
 - Goes an extra mile, you can start a webapp to manually curate your annotation.
 - Uses community driven databases and NCBI nr.

## Installation

You have multiple options to install and run decouphage:


Download and run standalone version:

    wget https://decouphage.com/decouphage
    chmod +x decouphage
    mv decouphage /usr/bin


Install with pip:

    pip install decouphage

Optional: Install dependencies in ubuntu:

    apt install ncbi-blast+ prodigal

Run with docker (Includes dependencies and databases):

    docker run decouphage/decouphage

## Databases

Decouphage database is derived from NCBI NR.


## Making custom databases

Make blast database

    makeblastdb -in database.fa -parse_seqids -blastdb_version 5 -dbtype prot


## Examples

Using blast and prodigal:

    $ decouphage --db phagedb ncbi --blast --prodigal --threads 8 genome.fasta

Using blast and phanotate:

    $ decouphage --db phagedb ncbi --scipy --phannotate --threads 8 genome.fasta

Using what is available and web curation:

    $ decouphage --db phagedb ncbi --curate --threads 8 genome.fasta

