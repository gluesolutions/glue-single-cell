from glue.viewers.matplotlib.state import (
    DeferredDrawCallbackProperty as DDCProperty,
    DeferredDrawSelectionCallbackProperty as DDSCProperty)
from glue.core.data_combo_helper import ComponentIDComboHelper, ComboHelper

from glue.viewers.scatter.state import ScatterViewerState

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

class QTLViewerState(ScatterViewerState):
    
    species = DDSCProperty(0, docstring='The species for displaying chromosome boundaries')
    pos_units = DDSCProperty(0, docstring='Units for gene and marker position')

    #chr_att = DDSCProperty(docstring='The attribute giving chromosome information')
    lod_att = DDSCProperty(docstring='The attribute giving the LOD score ', default_index=2)
    lod_thresh = DDCProperty(0, docstring='The LOD threshold for display and subsets')
    
    def __init__(self, **kwargs):
    
        super(QTLViewerState, self).__init__(**kwargs)
        
        #self.add_callback('layers', self._layers_changed)

        self.lod_att_helper = ComponentIDComboHelper(self, 'lod_att', numeric=True)
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

        self.pos_units_helper.choices = [1,1000,1_000_000,1_000_000_000]
        try:
            self.pos_units_helper.selection = self.pos_units_helper.choices[0]
        except IndexError:
            pass
        self.pos_units_helper.display = display_unit_names

        #self.chr_att_helper = ComponentIDComboHelper(self, 'chr_att', categorical=True, numeric=False) # Might be a problem if no X/Y/MT 

        #self.add_callback('lod_att', self._adjust_lod_thresh)

        #self.limits_cache = {}
        #                                            
        #self.cmap_lim_helper = StateAttributeLimitsHelper(self, attribute='lod_att',
        #                                                  lower='lod_vmin', upper='lod_vmax',
        #                                                  limits_cache=self.limits_cache)
        self.chr_pos = CHR_POSITIONS
        self.update_from_dict(kwargs)


    def _layers_changed(self, *args):
        
        layers_data = self.layers_data
        layers_data_cache = getattr(self, '_layers_data_cache', [])
        
        if layers_data == layers_data_cache:
            return
        
        self.x_att_helper.set_multiple_data(self.layers_data)
        self.y_att_helper.set_multiple_data(self.layers_data)
        self.lod_att_helper.set_multiple_data(self.layers_data)
        #self.chr_att_helper.set_multiple_data(self.layers_data)

        self._layers_data_cache = layers_data


    def _adjust_lod_thresh(self, *args):
        pass
        #print(f"{lod_att=}")

        # Not sure we need anything here... self.state.lod_thresh will be used by the viewer to apply_roi and by the layer_artist to show what to display, and that's it