import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from .constants import SAMPLES_COLOURS


def plot_position_graph(sample_variant_positions: dict[str, dict[str, list[int]]], ax: Axes | None = None) -> Axes:
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    # Build set of all unique chroms in sample variant positions
    all_chroms_set = set()
    for sample in sample_variant_positions:
        for chrom in sample_variant_positions[sample].keys():
            all_chroms_set.add(chrom)

    # Sort chromosomes
    all_chroms_sorted = sorted(all_chroms_set, key=lambda x: str(x).replace('chr', '').zfill(2))

    num_samples = len(sample_variant_positions)

    for i, chrom in enumerate(all_chroms_sorted):
        for j, sample in enumerate(sample_variant_positions):
            # Create 2 lists of coordinate positions.

            # X is the variants position in the chromosome
            x_positions = sample_variant_positions[sample][chrom]

            # Y is the position of the chromosome on the graph, offset by the sample
            y_position = ((i+1) * num_samples) - (j*.5)
            y_positions = [y_position] * len(x_positions)

            # Create label for first instance of this sample
            sample_label = sample if i == 0 else ""

            # Scatter plot
            ax.scatter(x_positions,
                       y_positions,
                       alpha=0.1,
                       s=5,
                       label=sample_label,
                       color=SAMPLES_COLOURS[j],
                       edgecolors='none'
                       )

    ax.set_yticks([i * num_samples for i in range(1, len(all_chroms_sorted)+1)], all_chroms_sorted)
    ax.set_xlabel("Genomic Position (bp)")
    ax.grid(axis='x', linestyle=':', alpha=0.6)
    ax.legend()

    return ax
