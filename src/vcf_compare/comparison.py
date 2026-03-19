import sys

from cyvcf2 import VCF
from matplotlib.axes import Axes
from simple_venn import venn4

from .models import VcfComparison


def filter_progress(fail_count, pass_count, last=False):
    end = '\r'
    # Print on every 10 pass records unless last
    if last or pass_count % 10 == 0:
        if last:
            end = '\n'
        print(f"Pass: {pass_count} Fail: {fail_count}", end=end, file=sys.stderr)

    return 1


def _vcf_records_to_sets(vcf: str) -> list[set[str]]:
    all = []
    passes = []

    pass_count = 0
    fail_count = 0

    for record in VCF(vcf):
        # TODO - this is where to add in selecting different INFO/FORMAT fields to compare
        record_str = f"{record.CHROM}\t{record.start}\t{record.end}\t{record.ALT[0]}"

        all.append(record_str)

        # FILTER is None if pass or . in cyvcf2
        if record.FILTER:
            fail_count += 1
        else:
            pass_count += 1
            passes.append(record_str)

            filter_progress(fail_count, pass_count)

    filter_progress(fail_count, pass_count, True)

    return [set(all), set(passes)]


def _original_venn_compare_sets(old_all, old_pass, new_all, new_pass):
    """  
    Original venn_compare subsets (15)
    [A, B, C, D, AB, AC, AD, BC, BD, CD, ABC, ABD, ACD, BCD, ABCD]
    """

    subsets = []
    
    subsets.append(  # Abcd
        old_all.difference(old_pass.union(new_all, new_pass))
    )
    subsets.append(  # aBcd
        old_pass.difference(old_all.union(new_all, new_pass))
    )
    subsets.append(  # abCd
        new_all.difference(old_all.union(old_pass, new_pass))
    )
    subsets.append(  # abcD
        new_pass.difference(old_all.union(old_pass, new_all))
    )
    subsets.append(  # ABcd
        old_all.intersection(old_pass).difference(new_all.union(new_pass))
    )
    subsets.append(  # AbCd
        old_all.intersection(new_all).difference(old_pass.union(new_pass))
    )
    subsets.append(  # AbcD
        old_all.intersection(new_pass).difference(old_pass.union(new_all))
    )
    subsets.append(  # aBCd
        old_pass.intersection(new_all).difference(old_all.union(new_pass))
    )
    subsets.append(  # aBcD
        old_pass.intersection(new_pass).difference(old_all.union(new_all))
    )
    subsets.append(  # abCD
        new_all.intersection(new_pass).difference(old_all.union(old_pass))
    )
    subsets.append(  # ABCd
        old_all.intersection(old_pass, new_all).difference(new_pass)
    )
    subsets.append(  # ABcD
        old_all.intersection(old_pass, new_pass).difference(new_all)
    )
    subsets.append(  # AbCD
        old_all.intersection(new_all, new_pass).difference(old_pass)
    )
    subsets.append(  # aBCD
        old_pass.intersection(new_all, new_pass).difference(old_all)
    )
    subsets.append(  # ABCD
        old_all.intersection(old_pass, new_all, new_pass)
    )

    return subsets


class Venn(VcfComparison):
    """ 
    Load 2 VCFs, plot unique/shared variants.
    TODO - how to configure which fields 
    TODO - setting for whether to look at PASS field
    """
    old_sets: list[set[str]]
    new_sets: list[set[str]]
    subsets: list


    def __init__(self, old_vcf: str, new_vcf: str) -> None:
        # Load VCFs
        print('Loading old: ' + old_vcf, file=sys.stderr)
        self.old_sets = _vcf_records_to_sets(old_vcf)

        print('Loading new: ' + new_vcf, file=sys.stderr)
        self.new_sets = _vcf_records_to_sets(new_vcf)
        
        # Generate subsets
        (old_all, old_pass, new_all, new_pass) = *self.old_sets, *self.new_sets

        print("Comparing sets...", file=sys.stderr)

        self.subsets = _original_venn_compare_sets(*self.old_sets, *self.new_sets)


    def plot(self) -> Axes:
        print("Plotting...", file=sys.stderr)
        subset_lens = []

        for e in self.subsets:
            subset_lens.append(len(e))

        ax = venn4(subsets=subset_lens,
                  set_labels=('old_a', 'old_p', 'new_a', 'new_p'),
                  set_label_fontsize=12,
                  subset_label_fontsize=10)
        
        ax.set_title("Venn of Old vs New")

        return ax


class Position(VcfComparison):
    """ Load multiple VCFs, plot each's variants positions () """
    def __init__(self) -> None:
        ...


class Metric(VcfComparison):
    """ Box plots of a specific metric - e.g. quality """
    def __init__(self) -> None:
        ...
