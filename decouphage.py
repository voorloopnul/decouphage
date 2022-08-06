#!/usr/bin/env python3
import json
import logging
import sys

import click
import time
from datetime import timedelta
from src.pipeline import Pipeline
from src.database import run_download

logger = logging.getLogger(__name__)


@click.command()
@click.option('--prodigal', is_flag=True, help='Use prodigal for orf calling instead of phanotate.')
@click.option('-d', '--db', 'database', type=click.Path(exists=True))
@click.option('-o', '--output', 'output', default='output.gbk')
@click.option('-t', '--threads', default=1, show_default=True)
@click.option('--tmpdir', 'tmp_dir', help="Folder for intermediate files.")
@click.option('--no_orf_calling', 'merge_gbk', is_flag=True, help="Annotate CDS from genbank file.")
@click.option('--locus_tag', 'locus_tag', default='dPHG', help="Locus tag prefix.")
@click.option('--download_db', 'download_db', is_flag=True, help="Download default database.")
@click.option('-v', '--verbose', is_flag=True, help="More verbose logging for debugging purpose.")
@click.argument('input_file', type=click.Path(exists=True))
def run_pipeline(prodigal, database, output, input_file, threads, tmp_dir, merge_gbk, locus_tag, download_db, verbose):
    if verbose:
        logging.basicConfig(level=logging.DEBUG)
        logger.setLevel(logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
        logger.setLevel(logging.INFO)

    if download_db:
        logger.error(" Option --download_db must be used alone.")
        sys.exit(1)

    logger.debug(json.dumps(locals(), indent=4, sort_keys=True))
    start_time = time.monotonic()

    Pipeline(database, input_file, prodigal, threads, output, tmp_dir, merge_gbk, locus_tag)
    end_time = time.monotonic()

    total_time = timedelta(seconds=end_time - start_time)
    logger.info(f"Phage annotation took {total_time}")


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == "--download_db":
        run_download()
    else:
        run_pipeline()
