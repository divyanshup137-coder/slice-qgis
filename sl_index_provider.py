"""Processing provider for the SLICE plugin."""

from qgis.core import QgsProcessingProvider
from .sl_index_algorithm import SLIndexAlgorithm


class SLIndexProvider(QgsProcessingProvider):
    """Registers the SLICE algorithm with the QGIS Processing framework."""

    def loadAlgorithms(self):
        self.addAlgorithm(SLIndexAlgorithm())

    def id(self):
        return "slice"

    def name(self):
        return "SLICE"

    def longName(self):
        return "SLICE: Stream Length-gradient Index by Constant Elevation"
