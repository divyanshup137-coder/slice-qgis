def classFactory(iface):
    from .sl_index_plugin import SLIndexPlugin
    return SLIndexPlugin(iface)
