import os

import click
from matplotlib import pyplot as plt

from . import __version__
from . comparison import Venn


"""py
parser = argparse.ArgumentParser(description='Compare 2 VCFs generating 4way venn.')
parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + getversion())
parser.add_argument('-o', '--old', dest='oldVcf', metavar='old.vcf[.gz]', help='VCF file for old result', required=True)
parser.add_argument('-n', '--new', dest='newVcf', metavar='new.vcf[.gz]', help='VCF file for new result', required=True)
parser.add_argument('-d', '--dir', dest='outDir', metavar='outDir', help='Directory for output', required=True)
parser.add_argument('-f', '--format', dest='imgFormat', metavar='[png|pdf|svg]', help='Format to save venn diagram', required=False, default='png')
parser.add_argument('-a', '--name', dest='name', metavar='sampleName', help='Prefix to name of output files', required=False, default='')
parser.add_argument('-y', '--info-filter', dest='infoFilters', metavar='INFO_TAG:MIN_PASS_VALUE:lt', action='append',
                    help='''Filter an INFO field entry. e.g. ASMD:140:lt treats positions with ASMD less than 140 as fail.
                            Valid comparators are lt (<), le (<=), gt (>), ge (>=), eq (==)''')
"""

@click.command()
@click.version_option(__version__)
@click.option("-o", "--out", type=click.Path(exists=True))
@click.option("-f", "--format", type=click.Choice(["png", "pdf", "svg"]), default="png")
@click.option("-i", "--info-filter", type=str)
@click.option("-a", "--name", type=str)
@click.argument("old", type=click.Path(exists=True, dir_okay=False))
@click.argument("new", type=click.Path(exists=True, dir_okay=False))
def venn_compare(out: str | None, format: str, info_filter: str, name: str, old: str, new: str):
    venn = Venn(old, new)
    
    figure = venn.plot()

    #plt.show()
    
    plt.plot()
    
    if out: 
        pre = f"{name}_" if name else ""
        plt.savefig(os.path.join(out, f"{pre}venn.{format}"))
