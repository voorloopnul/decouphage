#!/usr/bin/env python3
import json
import logging
import click
import time
from datetime import timedelta
from src.pipeline import Pipeline

logger = logging.getLogger(__name__)


@click.command()
@click.option('--prodigal', is_flag=True, help='Use prodigal for orf calling instead of phanotate.')
@click.option('-d', '--db', 'database', type=click.Path(exists=True))
@click.option('-o', '--output', 'output', default='output.gbk')
@click.option('-t', '--threads', default=1, show_default=True)
@click.option('--tmp_dir', 'tmp_dir', help="Folder for intermediate files.")
@click.option('--no-orf-calling', 'merge_gbk', is_flag=True, help="Annotate CDS from genbank file.")
@click.argument('input_file', type=click.Path(exists=True))
def run_pipeline(prodigal, database, output, input_file, threads, tmp_dir, merge_gbk):
    logger.info(json.dumps(locals(), indent=4, sort_keys=True))
    start_time = time.monotonic()

    Pipeline(database, input_file, prodigal, threads, output, tmp_dir, merge_gbk)
    end_time = time.monotonic()

    total_time = timedelta(seconds=end_time - start_time)
    logger.info(f"Phage annotation took {total_time}")


if __name__ == '__main__':
    run_pipeline()
