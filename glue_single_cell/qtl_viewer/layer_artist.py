from glue.viewers.scatter.layer_artist import ScatterLayerArtist
from glue.viewers.scatter.state import ScatterLayerState
from glue.viewers.scatter.python_export import python_export_scatter_layer
from glue.utils import defer_draw, broadcast_to, ensure_numerical
from glue.core.exceptions import IncompatibleAttribute

import numpy as np


class DensityMapLimits(object):

    contrast = 1
    
    def min(self, array):
        return 0
    
    def max(self, array):
        return 10. ** (np.log10(nanmax(array)) * self.contrast)


class QTLLayerArtist(ScatterLayerArtist):
    """
    This is not going to work for a density artist
    """
    
    _layer_state_cls = ScatterLayerState
    _python_exporter = python_export_scatter_layer
    
    def __init__(self, axes, viewer_state, layer_state=None, layer=None):
    
        super(QTLLayerArtist, self).__init__(axes, viewer_state,
                                                 layer_state=layer_state, layer=layer)
        
        # Watch for changes in the viewer state which would require the
        # layers to be redrawn
        self._viewer_state.add_global_callback(self._update_scatter)
        self.state.add_global_callback(self._update_scatter)
        
        # Scatter density
        self.density_auto_limits = DensityMapLimits()
        self._set_axes(axes)
        self.errorbar_index = 2
        self.vector_index = 3
        
        # NOTE: Matplotlib can't deal with NaN values in errorbar correctly, so
        # we need to prefilter values - the following variable is used to store
        # the mask for the values we keep, so that we can apply it to the color
        # See also https://github.com/matplotlib/matplotlib/issues/13799
        self._errorbar_keep = None

    @defer_draw
    def _update_data(self):
    
        print(f"{self._viewer_state.lod_att}")
        # Layer artist has been cleared already
        if len(self.mpl_artists) == 0:
            return
    
        try:
            if not self.state.density_map:
                x = ensure_numerical(self.layer[self._viewer_state.x_att].ravel())
    
        except (IncompatibleAttribute, IndexError):
            # The following includes a call to self.clear()
            self.disable_invalid_attributes(self._viewer_state.x_att)
            return
        else:
            self.enable()
    
        try:
            if not self.state.density_map:
                y = ensure_numerical(self.layer[self._viewer_state.y_att].ravel())
        except (IncompatibleAttribute, IndexError):
            # The following includes a call to self.clear()
            self.disable_invalid_attributes(self._viewer_state.y_att)
            return
        else:
            self.enable()
        
        if self._viewer_state.lod_att is not None and self._viewer_state.lod_thresh is not None:
            try:
                if not self.state.density_map:
                    lod = ensure_numerical(self.layer[self._viewer_state.lod_att].ravel())
            except (IncompatibleAttribute, IndexError):
                # The following includes a call to self.clear()
                self.disable_invalid_attributes(self._viewer_state.lod_att)
                return
            else:
                self.enable()
    
            lod_mask = np.ma.masked_where(lod < self._viewer_state.lod_thresh, lod)
            masked_x = np.ma.masked_where(np.ma.getmask(lod_mask), x)
            masked_y = np.ma.masked_where(np.ma.getmask(lod_mask), y)
        else:
            masked_x = x
            masked_y = y
    
        if self.state.markers_visible:
    
            if self.state.density_map:
                # We don't use x, y here because we actually make use of the
                # ability of the density artist to call a custom histogram
                # method which is defined on this class and does the data
                # access.
                self.plot_artist.set_data([], [])
                self.scatter_artist.set_offsets(np.zeros((0, 2)))
            else:
                self.density_artist.set_label(None)
                if self._use_plot_artist():
                    # In this case we use Matplotlib's plot function because it has much
                    # better performance than scatter.
                    self.plot_artist.set_data(masked_x, masked_y)
                else:
                    offsets = np.vstack((masked_x, masked_y)).transpose()
                    self.scatter_artist.set_offsets(offsets)
        else:
            self.plot_artist.set_data([], [])
            self.scatter_artist.set_offsets(np.zeros((0, 2)))