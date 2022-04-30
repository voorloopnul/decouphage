from src.pipeline import Pipeline
import click


@click.command()
@click.option('--prodigal', is_flag=True, help='Use prodigal for orf calling instead of phanotate.')
@click.option('-o', '--output', 'output', default='output.gbk')
@click.option('-t', '--threads', default=1, show_default=True)
@click.argument('input_file', type=click.Path(exists=True))
def run_pipeline(prodigal, output, input_file, threads):
    print(prodigal, output, input_file, threads)
    Pipeline(input_file, prodigal, threads, output)


if __name__ == '__main__':
    run_pipeline()
