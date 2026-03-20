from cyvcf2 import VCF
from matplotlib.axes import Axes
from simple_venn import venn2, venn4

from .models import VcfComparison
from .euler import plot_pass_fail_euler_diagram


def print_filter_progress(fail_count: int, pass_count: int, last: bool = False):
    """ Print pass and fail count on every 10 pass records, unless last """
    end = '\r'
    if last or pass_count % 10 == 0:
        if last:
            end = '\n'
        print(f"Pass: {pass_count} Fail: {fail_count}", end=end)
    return 1


def _vcf_records_to_pass_fail_sets(vcf_path: str, print_progress: bool = False) -> list[set[str]]:
    """ Read a VCF file and store information (TODO - what information?) about all and 
    pass-only/non-filtered variants, for comparison """
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


class VennVariantComparison(VcfComparison):
    """ Base class """
    old_sets: list[set[str]]
    new_sets: list[set[str]]
    sample_name: str

    def __init__(self, old_vcf_path: str, new_vcf_path: str, sample_name: str = "") -> None:
        # Load VCFs
        print('Loading old: ' + old_vcf_path)
        self.old_sets = _vcf_records_to_pass_fail_sets(old_vcf_path)

        print('Loading new: ' + new_vcf_path)
        self.new_sets = _vcf_records_to_pass_fail_sets(new_vcf_path)

        self.sample_name = sample_name


class Venn2VariantComparison(VennVariantComparison):
    """ Produces simple 2-way venn diagram of shared / unique variants """

    def __init__(self, old_vcf: str, new_vcf: str, sample_name: str = "", pass_only: bool = False):
        super().__init__(old_vcf, new_vcf, sample_name)
        self.pass_only = pass_only 

    def _2way_venn_compare_sets(self, old_set: set[str], new_set: set[str]) -> list[set[str]]:
        return [
            old_set.difference(new_set),    # A
            new_set.difference(old_set),    # B
            old_set.intersection(new_set)   # C
        ]

    def plot(self, ax: Axes | None = None) -> Axes:
        print("Comparing sets...")

        set_index = 1 if self.pass_only else 0
        subsets = self._2way_venn_compare_sets(self.old_sets[set_index], self.new_sets[set_index])

        print("Plotting...")
        subset_lens = []
        for subset in subsets:
            subset_lens.append(len(subset))

        ax = venn2(subsets=subset_lens,
                  set_labels=('old', 'new', 'shared'),
                  set_label_fontsize=12,
                  subset_label_fontsize=10,
                  ax=ax)
        
        prefix = "Venn of " if not self.sample_name else f"{self.sample_name} - "
        ax.set_title(prefix + "Old vs New")

        return ax


class Venn4VariantComparison(VennVariantComparison):
    """ Produces original 4-way vennCompare.py venn diagrams showing changes
    in pass variants """

    def _4way_venn_compare_sets(self,
                                    old_all: set[str],
                                    old_pass: set[str],
                                    new_all: set[str],
                                    new_pass: set[str]
                                ) -> list[set[str]]:
        """  
        Original venn_compare subsets (15)
        [A, B, C, D, AB, AC, AD, BC, BD, CD, ABC, ABD, ACD, BCD, ABCD]
        """

        return [
            old_all.difference(old_pass.union(new_all, new_pass)),              # Abcd
            old_pass.difference(old_all.union(new_all, new_pass)),              # aBcd
            new_all.difference(old_all.union(old_pass, new_pass)),              # abCd
            new_pass.difference(old_all.union(old_pass, new_all)),              # abcD
            old_all.intersection(old_pass).difference(new_all.union(new_pass)), # ABcd
            old_all.intersection(new_all).difference(old_pass.union(new_pass)), # AbCd
            old_all.intersection(new_pass).difference(old_pass.union(new_all)), # AbcD
            old_pass.intersection(new_all).difference(old_all.union(new_pass)), # aBCd
            old_pass.intersection(new_pass).difference(old_all.union(new_all)), # aBcD
            new_all.intersection(new_pass).difference(old_all.union(old_pass)), # abCD
            old_all.intersection(old_pass, new_all).difference(new_pass),       # ABCd
            old_all.intersection(old_pass, new_pass).difference(new_all),       # ABcD
            old_all.intersection(new_all, new_pass).difference(old_pass),       # AbCD
            old_pass.intersection(new_all, new_pass).difference(old_all),       # aBCD
            old_all.intersection(old_pass, new_all, new_pass)                   # ABCD
        ]

    def plot(self, ax: Axes | None = None) -> Axes:
        print("Comparing sets...")
        subsets = self._4way_venn_compare_sets(*self.old_sets, *self.new_sets)

        print("Plotting...")
        subset_lens = []
        for subset in subsets:
            subset_lens.append(len(subset))

        ax = venn4(subsets=subset_lens,
                   set_labels=('old_a', 'old_p', 'new_a', 'new_p'),
                   set_label_fontsize=12,
                   subset_label_fontsize=10,
                   ax=ax)
        
        prefix = "Venn of " if not self.sample_name else f"{self.sample_name} - "
        ax.set_title(prefix + "Old vs New")

        return ax


class EulerVariantComparison(VennVariantComparison):
    """ Produces pass / fail Euler diagram """

    def _euler_sets(self,
                    old_all: set[str],
                    old_pass: set[str],
                    new_all: set[str],
                    new_pass: set[str]
                    ) -> list[set[str]]:
        """ Generate 8 sets for euler diagram """
        return [
            old_all.difference(old_pass.union(new_all, new_pass)),          # 1 in old_all only
            new_all.difference(new_pass.union(old_all, old_pass)),          # 2 in new_all only 
            old_all.intersection(new_all).difference(old_pass, new_pass),   # 3 in old_all & new_all
            old_pass.difference(new_all),                                   # 4 in old_pass only
            new_pass.difference(old_all),                                   # 5 in new_pass only
            old_pass.intersection(new_all).difference(new_pass),            # 6 in old_pass & new_all (not new_pass)
            new_pass.intersection(old_all).difference(old_pass),            # 7 in new_pass & old_all (not old_pass)
            old_pass.intersection(new_pass)                                 # 8 in old_pass & new_pass
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
    """ Load multiple VCFs, plot each's variants positions () """
    def __init__(self) -> None:
        ...


class Metric(VcfComparison):
    """ Box plots of a specific metric - e.g. quality """
    def __init__(self) -> None:
        ...
