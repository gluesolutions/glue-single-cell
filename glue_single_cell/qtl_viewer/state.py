from glue.viewers.matplotlib.state import (
    DeferredDrawCallbackProperty as DDCProperty,
    DeferredDrawSelectionCallbackProperty as DDSCProperty)
from glue.core.data_combo_helper import ComponentIDComboHelper

from glue.viewers.scatter.state import ScatterViewerState

__all__ = ['QTLViewerState']

class QTLViewerState(ScatterViewerState):
    
    lod_att = DDSCProperty(docstring='The attribute giving the LOD score ', default_index=2)
    lod_thresh = DDCProperty(0, docstring='The LOD threshold for display and subsets')
    
    def __init__(self, **kwargs):
    
        super(QTLViewerState, self).__init__(**kwargs)
        
        #self.add_callback('layers', self._layers_changed)

        self.lod_att_helper = ComponentIDComboHelper(self, 'lod_att', numeric=True)
        #self.add_callback('lod_att', self._adjust_lod_thresh)

        #self.limits_cache = {}
        #                                            
        #self.cmap_lim_helper = StateAttributeLimitsHelper(self, attribute='lod_att',
        #                                                  lower='lod_vmin', upper='lod_vmax',
        #                                                  limits_cache=self.limits_cache)
        self.update_from_dict(kwargs)


    def _layers_changed(self, *args):
        
        layers_data = self.layers_data
        layers_data_cache = getattr(self, '_layers_data_cache', [])
        
        if layers_data == layers_data_cache:
            return
        
        self.x_att_helper.set_multiple_data(self.layers_data)
        self.y_att_helper.set_multiple_data(self.layers_data)
        self.lod_att_helper.set_multiple_data(self.layers_data)

        self._layers_data_cache = layers_data


    def _adjust_lod_thresh(self, *args):
        pass
        #print(f"{lod_att=}")

        # Not sure we need anything here... self.state.lod_thresh will be used by the viewer to apply_roi and by the layer_artist to show what to display, and that's it