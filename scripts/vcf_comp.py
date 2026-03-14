import click

from vcf_compare import __version__


@click.command()
@click.version_option(__version__)
def vcf_compare():
    """ Command-line tool for VCF comparison """
    raise NotImplementedError


if __name__ == "__main__":
    vcf_compare()
