# Decouphage: the art of decorating a Phage genome by gluing colored gene labels onto it.

## Description

 - **Phages only** as the name suggests decouphage is a tool designed to annotate phage genomes. This decision has been made
in order to tweak the annotation process regarding tools, preferences & databases and to streamline further development &
maintenance of the software.


## Installation

Install with pip:

    pip install decouphage

Install dependencies in ubuntu:

    apt install ncbi-blast+ prodigal

## Making custom databases

Make blast database

    makeblastdb -in database.fa -parse_seqids -blastdb_version 5 -dbtype prot


## Examples

Using blast and prodigal:

    $ decouphage --db phagedb ncbi --align blast --orf prodigal --threads 8 genome.fasta

Using blast and phanotate:

    $ decouphage --db phagedb ncbi --align scipy --orf phannotate --threads 8 genome.fasta
