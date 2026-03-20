import os

import click
from matplotlib import pyplot as plt

from . import __version__
from . comparison import EulerVariantComparison


@click.command()
@click.version_option(__version__)
@click.option("-o", "--out", type=click.Path(exists=True))
@click.option("-f", "--format", type=click.Choice(["png", "pdf", "svg"]), default="png")
@click.option("-i", "--info-filter", type=str)
@click.option("-a", "--name", type=str)
@click.argument("old", type=click.Path(exists=True, dir_okay=False))
@click.argument("new", type=click.Path(exists=True, dir_okay=False))
def venn_compare(out: str | None, format: str, info_filter: str, name: str, old: str, new: str):
    venn = EulerVariantComparison(old, new, name)
    
    ax = venn.plot()

    # TODO - is this needed for CLI mode?
    plt.plot()
    
    if out: 
        pre = f"{name}_" if name else ""
        plt.savefig(os.path.join(out, f"{pre}venn.{format}"))
