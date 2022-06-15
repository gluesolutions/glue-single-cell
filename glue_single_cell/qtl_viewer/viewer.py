from glue.viewers.scatter.qt.data_viewer import ScatterViewer
from glue.viewers.scatter.qt.layer_style_editor import ScatterLayerStyleEditor
from glue.viewers.scatter.layer_artist import ScatterLayerArtist
from glue.utils import defer_draw, decorate_all_methods

from .layer_artist import QTLLayerArtist
from .qt.options_widget import QTLOptionsWidget
from .state import QTLViewerState

__all__ = ['QTLViewer']

@decorate_all_methods(defer_draw)
class QTLViewer(ScatterViewer):
    LABEL = 'QTL Viewer'
    _layer_style_widget_cls = ScatterLayerStyleEditor # We can just reuse this layer style for now, although eventually we should trim options that do not make sense
    _state_cls = QTLViewerState
    _options_cls = QTLOptionsWidget # This has the LOD level widget
    _data_artist_cls = QTLLayerArtist 
    _subset_artist_cls = QTLLayerArtist
    
    tools = ['select:rectangle', 'select:xrange',
     'select:yrange']

    def __init__(self, session, parent=None, state=None):
        super(QTLViewer, self).__init__(session, parent=parent, state=state)
