"""QGIS plugin entry point: registers and unregisters the SLICE provider."""

from qgis.core import QgsApplication
from .sl_index_provider import SLIndexProvider


class SLIndexPlugin:
    """Adds the SLICE processing provider on load and removes it on unload."""

    def __init__(self, iface):
        self.iface = iface
        self.provider = None

    def initGui(self):
        self.provider = SLIndexProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

    def unload(self):
        QgsApplication.processingRegistry().removeProvider(self.provider)
