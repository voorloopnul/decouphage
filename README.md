![Decouphage logo](https://raw.githubusercontent.com/voorloopnul/voorloopnul/357a7ead62584e352c61b008790fe38d4aff5664/logos/decouphage.png)


# Decouphage: the art of decorating a Phage genome by gluing feature cutouts into it.

### Relevant branches

 - [main branch](https://github.com/voorloopnul/decouphage/tree/main): stable version available in [pypi](https://pypi.org/project/decouphage/) and [dockerhub](https://hub.docker.com/r/voorloop/decouphage).
 - [dev branch](https://github.com/voorloopnul/decouphage/tree/dev): development branch with most recent version of decouphage.

## Table of contents

- [Decouphage: the art of decorating a Phage genome by gluing feature cutouts into it.](#decouphage--the-art-of-decorating-a-phage-genome-by-gluing-feature-cutouts-into-it)
    + [Relevant branches](#relevant-branches)
  * [Description](#description)
    + [Highlights](#highlights)
  * [Methods and Validation](#methods-and-validation)
    + [Validation](#validation)
  * [How can I use decouphage](#how-can-i-use-decouphage)
    + [Options](#options)
    + [I want to discover and annotate a lot of ORFs](#i-want-to-discover-and-annotate-a-lot-of-orfs)
    + [I want to use prodigal to find my genes](#i-want-to-use-prodigal-to-find-my-genes)
    + [I have a genbank with poor annotation and want more](#i-have-a-genbank-with-poor-annotation-and-want-more)
  * [Installation](#installation)
    + [Ubuntu](#ubuntu)
    + [Docker](#docker)
  * [Databases](#databases)
    + [Downloading database](#downloading-database)
    + [Making custom databases](#making-custom-databases)

## Description

As the name suggests decouphage is a tool designed to annotate phage genomes. It only external dependency is ncbi-blast+
everything else is optional. 
 
### Highlights

 - Can be easily installed in Linux or Mac computers. Only requirement is ncbi-blast+.
 - Can be extended with prodigal, but as default it uses phanotate for ORF calling. 
 - Decouphage is fast, using a Macbook **most phage genomes can be annotated in less than a minute**.
 - Uses ncbi NR database containing non-identical sequences from GenBank CDS translations, PDB, Swiss-Prot, PIR, and PRF. 
 

## Methods and Validation

### Methods

To be written...

### Validation

Decouphage validation was made in comparison to RAST(Rapid Annotation using Subsystem Technology), a tool that is often
praised for its good Prokaryotic annotation capabilities.

Decouphage outperforms RAST when calling some of the most relevant product categories:

![alt text](validation/decouphage_image_01.png?raw=true)

The CDS annotation agreement between Decouphage and RAST is high, reaching up to 94% for some products:

| Enzyme | Agreement rate with RAST |
| ------ | ------------------------- |
| endonuclease |  94% |
| exonuclease | 58% |
| helicase | 70% |
| hydrolase | 73% |
| kinase |  86% |
| ligase |  94% |
| methyltransferase | 65% |
| polymerase | 76% |
| primase | 78% |
| protease | 85% |
| recombinase | 28% |
| reductase | 90% |
| synthase | 84% |
| terminase | 94% |
| transferase | 60% |

A precise comparison of product-to-position is difficult given differences in spelling, typos, synonyms, and 
interchangeable names, but the table above can give a good sense of the similarities.

To corroborate the surplus of annotations that decouphage achieves, the amount of "hypothetical protein" and 
"Phage protein" was also checked:

| Product | Decouphage | Rast | Agreement rate with RAST |
| ------- | ---------- | ---- | ------------------------ |
| hypothetical protein | 3945 | 6302 | 53% |
| phage protein |    0 | 1626 | N/A<sup>1</sup> |
| Total products |    9692 | 9692 | N/A<sup>2</sup>  |

1. Decouphage does not include products containing "phage protein" as they usually are a noise source.
2. The genbank file generated by RAST was used as input for decouphage to ensure no difference in the number of CDS.

This table shows that Decouphage potentially assigns 2x more meaningful products than RAST when annotating a phage 
genome.

**DISCLAIMER: While decouphage validation appears to outperform RAST when annotating phages, remember that RAST is - 
by far - a more advanced tool with consolidated methods. The results shown here are more a reflection of the database 
used than the underlying algorithms.** 


## How can I use decouphage

### Options

    Usage: decouphage [OPTIONS] INPUT_FILE
    
    Options:
      --prodigal             Use prodigal for orf calling instead of phanotate.
      -d, --db PATH
      -o, --output TEXT
      -t, --threads INTEGER  [default: 1]
      --tmpdir TEXT          Folder for intermediate files.
      --no_orf_calling       Annotate CDS from genbank file.
      --locus_tag TEXT       Locus tag prefix.
      --download_db          Download default database.
      -v, --verbose          More verbose logging for debugging purpose.
      --help                 Show this message and exit.


### I want to discover and annotate a lot of ORFs
 
    decouphage genome.fasta -o genome.gb

### I want to use prodigal to find my genes

    decouphage genome.fasta -o genome.gb --prodigal

### I have a genbank with poor annotation and want more

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

### Making custom databases

Make blast database

    makeblastdb -in database.fa -parse_seqids -blastdb_version 5 -dbtype prot
