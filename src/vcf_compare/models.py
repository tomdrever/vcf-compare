from abc import ABC

from matplotlib.axes import Axes

class VcfComparison(ABC):
    def plot(self) -> Axes:
        """ Base function to produce plot to show """
        ...
