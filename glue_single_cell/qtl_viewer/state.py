import numpy as np

from glue.core import BaseData, Subset

from glue.viewers.matplotlib.state import (MatplotlibDataViewerState,
                                           DeferredDrawCallbackProperty as DDCProperty,
                                           DeferredDrawSelectionCallbackProperty as DDSCProperty)
from glue.core.data_combo_helper import ComponentIDComboHelper, ComboHelper
from glue.core.state_objects import StateAttributeLimitsHelper

#from glue.viewers.scatter.state import ScatterViewerState

__all__ = ['QTLViewerState']

CHR_POSITIONS = {
    'Mouse':{ # http://www.informatics.jax.org/mgihome/other/mouse_facts1.shtml
        'Names':['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X','Y','MT'],
        'Lengths':[195,182,160,157,152,150,145,130,124,131,122,120,121,125,104,98,95,91,61,169,91,0.01],
        'GridSize':200_000_000, # Genome position seems to use a fixed grid so that there is emtpy space in later chromosomes
    },
    'Human':{ # http://www.insilicase.com/Web/Chromlen.aspx
        'Names':['1','2','3','4','5','6','7','8','9','10',
                    '11','12','13','14','15','16','17','18','19','20','21','22','22','X','Y','MT'],
        'GridSize':250_000_000, #This is a guess
                    
    }
    
}
UNITS_LOOKUP = {1:'bp',1000:'kb',1_000_000:'Mb',1_000_000_000:'Gb'}

class QTLViewerState(MatplotlibDataViewerState):

    x_att = DDSCProperty(docstring='The attribute to show on the x-axis', default_index=0)
    y_att = DDSCProperty(docstring='The attribute to show on the y-axis', default_index=1)
    dpi = DDCProperty(72, docstring='The resolution (in dots per inch) of density maps, if present')

    species = DDSCProperty(0, docstring='The species for displaying chromosome boundaries')
    pos_units = DDSCProperty(0, docstring='Units for gene and marker position')

    lod_att = DDSCProperty(docstring='The attribute giving the LOD score ', default_index=2)
    lod_thresh = DDCProperty(0, docstring='The LOD threshold for display and subsets')
    
    plot_mode = DDSCProperty(docstring="Whether to plot the data in cartesian, polar or another projection")
    angle_unit = DDSCProperty(docstring="Whether to use radians or degrees for any angular coordinates")

    def __init__(self, **kwargs):
        super().__init__()

        self.limits_cache = {}

        self.x_lim_helper = StateAttributeLimitsHelper(self, attribute='x_att',
                                                       lower='x_min', upper='x_max',
                                                       log='x_log', margin=0.04,
                                                       limits_cache=self.limits_cache)

        self.y_lim_helper = StateAttributeLimitsHelper(self, attribute='y_att',
                                                       lower='y_min', upper='y_max',
                                                       log='y_log', margin=0.04,
                                                       limits_cache=self.limits_cache)

        self.add_callback('layers', self._layers_changed)
        
        self.x_att_helper = ComponentIDComboHelper(self, 'x_att', pixel_coord=True, world_coord=True)
        self.y_att_helper = ComponentIDComboHelper(self, 'y_att', pixel_coord=True, world_coord=True)
        self.lod_att_helper = ComponentIDComboHelper(self, 'lod_att', numeric=True, categorical=False)
    
        self.plot_mode_helper = ComboHelper(self, 'plot_mode')
        self.plot_mode_helper.choices = ['rectilinear']
        self.plot_mode_helper.selection = 'rectilinear'

        self.angle_unit_helper = ComboHelper(self, 'angle_unit')
        self.angle_unit_helper.choices = ['radians', 'degrees']
        self.angle_unit_helper.selection = 'radians'

        self.species_helper = ComboHelper(self, 'species')
        self.pos_units_helper = ComboHelper(self, 'pos_units')

        self.chr_pos = CHR_POSITIONS

        self.species_helper.choices = list(self.chr_pos.keys())
        try:
            self.species_helper.selection = self.species_helper.choices[0]
        except IndexError:
            pass
        self.species_helper.display = str


        def display_unit_names(unit):
            return UNITS_LOOKUP[unit]

        self.pos_units_helper.choices = [1, 1000, 1_000_000, 1_000_000_000]
        try:
            self.pos_units_helper.selection = self.pos_units_helper.choices[0]
        except IndexError:
            pass
        self.pos_units_helper.display = display_unit_names

        self.chr_pos = CHR_POSITIONS

        self.update_from_dict(kwargs)

        self.add_callback('x_log', self._reset_x_limits)
        self.add_callback('y_log', self._reset_y_limits)


    def _layers_changed(self, *args):
        
        layers_data = self.layers_data
        layers_data_cache = getattr(self, '_layers_data_cache', [])
        
        if layers_data == layers_data_cache:
            return
        
        self.x_att_helper.set_multiple_data(self.layers_data)
        self.y_att_helper.set_multiple_data(self.layers_data)
        self.lod_att_helper.set_multiple_data(self.layers_data)

        self._layers_data_cache = layers_data

    # Everything below here is copy-paste from ScatterViewerState

    def _reset_x_limits(self, *args):
        if self.x_att is None:
            return
        self.x_lim_helper.percentile = 100
        self.x_lim_helper.update_values(force=True)

    def _reset_y_limits(self, *args):
        if self.y_att is None:
            return
        self.y_lim_helper.percentile = 100
        self.y_lim_helper.update_values(force=True)

    def reset_limits(self):
        if not self.using_polar:
            self._reset_x_limits()
        self._reset_y_limits()

    def flip_x(self):
        self.x_lim_helper.flip_limits()

    def flip_y(self):
        self.y_lim_helper.flip_limits()

    @property
    def x_categories(self):
        return self._categories(self.x_att)

    @property
    def y_categories(self):
        return self._categories(self.y_att)

    def _categories(self, cid):

        categories = []

        for layer_state in self.layers:

            if isinstance(layer_state.layer, BaseData):
                layer = layer_state.layer
            else:
                layer = layer_state.layer.data

            try:
                if layer.data.get_kind(cid) == 'categorical':
                    categories.append(layer.data.get_data(cid).categories)
            except IncompatibleAttribute:
                pass

        if len(categories) == 0:
            return None
        else:
            return np.unique(np.hstack(categories))

    @property
    def x_kinds(self):
        return self._component_kinds(self.x_att)

    @property
    def y_kinds(self):
        return self._component_kinds(self.y_att)

    def _component_kinds(self, cid):

        # Construct list of component kinds over all layers

        kinds = set()

        for layer_state in self.layers:

            if isinstance(layer_state.layer, BaseData):
                layer = layer_state.layer
            else:
                layer = layer_state.layer.data

            try:
                kinds.add(layer.data.get_kind(cid))
            except IncompatibleAttribute:
                pass

        return kinds

    def _adjust_lod_thresh(self, *args):
        pass
        #print(f"{lod_att=}")

        # Not sure we need anything here... self.state.lod_thresh will be used by the viewer to apply_roi and by the layer_artist to show what to display, and that's it