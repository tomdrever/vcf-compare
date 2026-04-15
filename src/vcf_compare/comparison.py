from typing import Any, Callable
from collections import defaultdict

from cyvcf2 import VCF, Variant
from matplotlib.axes import Axes
from simple_venn import venn2, venn4
import matplotlib.pyplot as plt

from .models import VcfComparison
from .euler import plot_pass_fail_euler_diagram
from .position import plot_position_graph


def print_filter_progress(fail_count: int, pass_count: int, last: bool = False):
    """Print pass and fail count on every 10 pass records, unless last"""
    end = "\r"
    if last or pass_count % 10 == 0:
        if last:
            end = "\n"
        print(f"Pass: {pass_count} Fail: {fail_count}", end=end)
    return 1


def _vcf_records_to_pass_fail_sets(vcf_path: str, print_progress: bool = False) -> list[set[str]]:
    """Read a VCF file and store information (TODO - what information?) about all and
    pass-only/non-filtered variants, for comparison"""
    all = []
    passes = []

    pass_count = 0
    fail_count = 0

    for record in VCF(vcf_path):
        # TODO - this is where to add in selecting different INFO/FORMAT fields to compare
        record_str = f"{record.CHROM}\t{record.start}\t{record.end}\t{record.ALT[0]}"

        all.append(record_str)

        # FILTER is None if pass or . in cyvcf2
        if record.FILTER:
            fail_count += 1
        else:
            pass_count += 1
            passes.append(record_str)

            if print_progress:
                print_filter_progress(fail_count, pass_count)

    if print_progress:
        print_filter_progress(fail_count, pass_count, True)

    return [set(all), set(passes)]


def _get_info_field_metric(record: Variant, field: str) -> int | float:
    metric = record.INFO[field]
    if type(metric) not in [int, float]:
        raise ValueError(f"Error parsing metric INFO field '{field}' - invalid type {type(metric)}")
    return metric


def _resolve_metric(metric: str | Callable[[Any], Any]) -> Callable[[Any], Any]:
    if callable(metric):
        return metric
    if metric == "QUAL":
        return lambda r: r.QUAL
    if metric.startswith("INFO."):
        field = metric[5:]
        return lambda r: _get_info_field_metric(r, field)
    raise ValueError(f"Unknown metric {metric!r}. Use 'QUAL', 'INFO.<field>', or a callable.")


def _vcf_records_to_metric_list(
    vcf_path: str,
    pass_only: bool,
    metric: str | Callable[[Any], Any]
) -> list[Any]:
    metrics = []
    extractor = _resolve_metric(metric)

    for record in VCF(vcf_path):
        # If pass_only = true only add if no record.FILTER
        if not pass_only or not record.FILTER:
            metrics.append(extractor(record))

    return metrics


class VennVariantComparison(VcfComparison):
    """Base class for Old vs New comparisons"""
    old_sets: list[set[str]]
    new_sets: list[set[str]]
    sample_name: str
    a_name: str
    b_name: str

    def __init__(
        self,
        old_vcf_path: str,
        new_vcf_path: str,
        sample_name: str = "",
        a_name: str = "Old",
        b_name: str = "New",
    ) -> None:
        # Load VCFs
        print("Loading old: " + old_vcf_path)
        self.old_sets = _vcf_records_to_pass_fail_sets(old_vcf_path)

        print("Loading new: " + new_vcf_path)
        self.new_sets = _vcf_records_to_pass_fail_sets(new_vcf_path)

        self.sample_name = sample_name
        self.a_name = a_name
        self.b_name = b_name


class Venn2VariantComparison(VennVariantComparison):
    """Produces simple 2-way venn diagram of shared / unique variants"""

    def __init__(
        self,
        old_vcf: str,
        new_vcf: str,
        sample_name: str = "",
        a_name: str = "Old",
        b_name: str = "New",
        pass_only: bool = False,
    ) -> None:
        super().__init__(old_vcf, new_vcf, sample_name, a_name, b_name)
        self.pass_only = pass_only

    def _2way_venn_compare_sets(self, old_set: set[str], new_set: set[str]) -> list[set[str]]:
        return [
            old_set.difference(new_set),  # 1. in old but not new
            new_set.difference(old_set),  # 2. in new but not old
            old_set.intersection(new_set),  # 3. in both old and new
        ]

    def plot(self, ax: Axes | None = None, title: str | None = None) -> Axes:
        print("Comparing sets...")

        set_index = 1 if self.pass_only else 0
        subsets = self._2way_venn_compare_sets(self.old_sets[set_index], self.new_sets[set_index])

        print("Plotting...")
        subset_lens = []
        for subset in subsets:
            subset_lens.append(len(subset))

        ax = venn2(
            subsets=subset_lens,
            set_labels=(self.a_name, self.b_name, "shared"),
            set_label_fontsize=18,
            subset_label_fontsize=16,
            set_colors=["#7B6FCF", "#2A9D7A"],
            ax=ax,
        )

        prefix = "Venn of " if not self.sample_name else f"{self.sample_name} - "

        if not title:
            ax.set_title(prefix + f"{self.a_name} vs {self.b_name}", fontdict={"fontsize": 18})
        else:
            ax.set_title(title, fontdict={"fontsize": 18})

        return ax


class Venn4VariantComparison(VennVariantComparison):
    """Produces original 4-way vennCompare.py venn diagrams showing changes
    in pass variants"""

    def _4way_venn_compare_sets(
        self, old_all: set[str], old_pass: set[str], new_all: set[str], new_pass: set[str]
    ) -> list[set[str]]:
        """
        Original venn_compare subsets (15)
        [A, B, C, D, AB, AC, AD, BC, BD, CD, ABC, ABD, ACD, BCD, ABCD]
        """

        return [
            old_all.difference(old_pass.union(new_all, new_pass)),  # Abcd
            old_pass.difference(old_all.union(new_all, new_pass)),  # aBcd
            new_all.difference(old_all.union(old_pass, new_pass)),  # abCd
            new_pass.difference(old_all.union(old_pass, new_all)),  # abcD
            old_all.intersection(old_pass).difference(new_all.union(new_pass)),  # ABcd
            old_all.intersection(new_all).difference(old_pass.union(new_pass)),  # AbCd
            old_all.intersection(new_pass).difference(old_pass.union(new_all)),  # AbcD
            old_pass.intersection(new_all).difference(old_all.union(new_pass)),  # aBCd
            old_pass.intersection(new_pass).difference(old_all.union(new_all)),  # aBcD
            new_all.intersection(new_pass).difference(old_all.union(old_pass)),  # abCD
            old_all.intersection(old_pass, new_all).difference(new_pass),  # ABCd
            old_all.intersection(old_pass, new_pass).difference(new_all),  # ABcD
            old_all.intersection(new_all, new_pass).difference(old_pass),  # AbCD
            old_pass.intersection(new_all, new_pass).difference(old_all),  # aBCD
            old_all.intersection(old_pass, new_all, new_pass),  # ABCD
        ]

    def plot(self, ax: Axes | None = None) -> Axes:
        print("Comparing sets...")
        subsets = self._4way_venn_compare_sets(*self.old_sets, *self.new_sets)

        print("Plotting...")
        subset_lens = []
        for subset in subsets:
            subset_lens.append(len(subset))

        ax = venn4(
            subsets=subset_lens,
            set_labels=("old_a", "old_p", "new_a", "new_p"),
            set_label_fontsize=12,
            subset_label_fontsize=10,
            ax=ax,
        )

        prefix = "Venn of " if not self.sample_name else f"{self.sample_name} - "
        ax.set_title(prefix + "Old vs New")

        return ax


class EulerVariantComparison(VennVariantComparison):
    """Produces pass / fail Euler diagram"""

    def _euler_sets(
        self, old_all: set[str], old_pass: set[str], new_all: set[str], new_pass: set[str]
    ) -> list[set[str]]:
        """Generate 8 sets for euler diagram"""
        return [
            old_all.difference(old_pass.union(new_all, new_pass)),  # 1 in old_all only
            new_all.difference(new_pass.union(old_all, old_pass)),  # 2 in new_all only
            old_all.intersection(new_all).difference(old_pass, new_pass),  # 3 in old_all & new_all
            old_pass.difference(new_all),  # 4 in old_pass only
            new_pass.difference(old_all),  # 5 in new_pass only
            old_pass.intersection(new_all).difference(new_pass),  # 6 in old_pass & new_all (not new_pass)
            new_pass.intersection(old_all).difference(old_pass),  # 7 in new_pass & old_all (not old_pass)
            old_pass.intersection(new_pass),  # 8 in old_pass & new_pass
        ]

    def plot(self, ax: Axes | None = None) -> Axes:
        print("Comparing sets...")
        subsets = self._euler_sets(*self.old_sets, *self.new_sets)

        print("Plotting...")
        subset_lens = []
        for subset in subsets:
            subset_lens.append(len(subset))

        ax = plot_pass_fail_euler_diagram(subset_lens, ax=ax)

        prefix = "Venn of " if not self.sample_name else f"{self.sample_name} - "
        ax.set_title(prefix + "Old vs New")

        return ax


class Position(VcfComparison):
    """Load multiple VCFs, plot each's variants positions"""

    variant_positions: dict[str, set[str]]

    def __init__(
        self,
        *vcfs: tuple[str, str],
        pass_only: bool = False,
        unique_only: bool = False
    ) -> None:
        """Supports multiple VCFs via args, which are all assumed to be tuples in format [name, path]"""

        self.variant_positions = {}

        # Load VCFs
        for vcf in vcfs:
            name = vcf[0]
            path = vcf[1]

            vcf_chroms = set()
            for record in VCF(path):
                # If pass-only, filter
                if not pass_only or not record.FILTER:
                    vcf_chroms.add(f"{record.CHROM}-{record.POS}")

            self.variant_positions[name] = vcf_chroms

        if unique_only:
            # Update variant positions to be only those not found in other variants
            for sample in self.variant_positions:
                sample_variants = self.variant_positions[sample]

                # Get all variants from other samples
                other_variants = set()
                for other_sample in self.variant_positions:
                    if sample != other_sample:
                        other_variants.update(self.variant_positions[other_sample])

                self.variant_positions[sample] = sample_variants.difference(other_variants)

    def plot(self, ax: Axes | None = None) -> Axes:
        # Split variant position by chromosome
        variant_positions_by_chr: dict[str, dict[str, list[int]]] = {}

        for sample in self.variant_positions:
            sample_var_pos_by_chr = defaultdict(list)
            for variant_position in self.variant_positions[sample]:
                chr, position = variant_position.split("-")
                sample_var_pos_by_chr[chr].append(int(position))

            variant_positions_by_chr[sample] = sample_var_pos_by_chr

        return plot_position_graph(variant_positions_by_chr, ax)


class Metric(VcfComparison):
    """ Box plots of a specific metric for multiple VCFs - e.g. quality """

    metrics: dict[str, list[Any]] # dict of {sample : variant_metric_values}

    def __init__(
        self,
        *vcfs: tuple[str, str],
        pass_only: bool = False,
        metric: str | Callable[[Any], Any],
    ) -> None:
        """Supports multiple VCFs via args, which are all assumed to be tuples in format [name, path]"""
        self.metrics = {}
        self.metric = metric
        self.pass_only = pass_only

        # Load VCFs
        for vcf in vcfs:
            name = vcf[0]
            path = vcf[1]

            print(f"Loading {name}: {path}")
            self.metrics[name] = _vcf_records_to_metric_list(path, pass_only, metric)


    def plot(self, ax: Axes | None = None) -> Axes:
        if ax is None:
            _, ax = plt.subplots(figsize=(10, 6))

        metrics_lists = [self.metrics[sample] for sample in self.metrics]

        ax.boxplot(metrics_lists, labels=self.metrics.keys(), patch_artist=True)
        ax.set_xticklabels(list(self.metrics.keys()))

        metric_label = self.metric if isinstance(self.metric, str) else "custom metric"
        pass_label = "(passes only) " if self.pass_only else ""
        ax.set_xlabel("Variant set")
        ax.set_ylabel(metric_label)
        ax.set_title(f"{metric_label} {pass_label} by sample")
        ax.legend()

        return ax
