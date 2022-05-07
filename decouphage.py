import json
import logging
import click
import time
from datetime import timedelta
from src.pipeline import Pipeline

logger = logging.getLogger(__name__)


@click.command()
@click.option('--prodigal', is_flag=True, help='Use prodigal for orf calling instead of phanotate.')
@click.option('-o', '--output', 'output', default='output.gbk')
@click.option('-t', '--threads', default=1, show_default=True)
@click.argument('input_file', type=click.Path(exists=True))
def run_pipeline(prodigal, output, input_file, threads):
    logger.info(json.dumps(locals(), indent=4, sort_keys=True))
    start_time = time.monotonic()

    Pipeline(input_file, prodigal, threads, output)
    end_time = time.monotonic()

    total_time = timedelta(seconds=end_time - start_time)
    logger.info(f"Phage annotaion took {total_time}")


if __name__ == '__main__':
    run_pipeline()
