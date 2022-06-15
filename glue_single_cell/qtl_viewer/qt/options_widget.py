import os
from qtpy import QtWidgets
from echo.qt import autoconnect_callbacks_to_qt
from glue.utils.qt import load_ui, fix_tab_widget_fontsize


__all__ = ['QTLOptionsWidget']


class QTLOptionsWidget(QtWidgets.QWidget):

    def __init__(self, viewer_state, session, parent=None):
    
        super(QTLOptionsWidget, self).__init__(parent=parent)
    
        self.ui = load_ui('options_widget.ui', self,
                          directory=os.path.dirname(__file__))
    
        fix_tab_widget_fontsize(self.ui.tab_widget)
        connect_kwargs = {'lod_thresh': dict(value_range=(1,200))}  # We might want to take these from data
    
        self._connections = autoconnect_callbacks_to_qt(viewer_state, self.ui, connect_kwargs)
        connect_kwargs = {'alpha': dict(value_range=(0, 1))}
        self._connections_legend = autoconnect_callbacks_to_qt(viewer_state.legend, self.ui.legend_editor.ui, connect_kwargs)
    
        self.viewer_state = viewer_state
