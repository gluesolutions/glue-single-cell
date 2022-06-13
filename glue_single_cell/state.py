from glue.core.state_objects import State
from echo import SelectionCallbackProperty, CallbackProperty
from glue.core.data_combo_helper import DataCollectionComboHelper, ComponentIDComboHelper, ComboHelper, ManualDataComboHelper

__all__ = ['DiffGeneExpState', 'PCASubsetState', 'GSEApyState']


class GSEApyState(State):
    data = SelectionCallbackProperty() # The data object to enrich
    subset = SelectionCallbackProperty() # Optionally a subset that can be applied to data
    #organism = SelectionCallbackProperty()
    gene_set = SelectionCallbackProperty() # The Enrichr library to use. Could be linked to organism to get a filtered list.
    # This is a terrible name, because it is confusing with gene subsets, but this is the GSEApy parameter name
    gene_att = SelectionCallbackProperty() # The attribute with gene labels in it

    
    def __init__(self, data_collection):
        
        super(GSEApyState, self).__init__()

        self.data_collection = data_collection
        self.data_helper = DataCollectionComboHelper(self, 'data', data_collection)
        self.subset_helper = ComboHelper(self, 'subset')
        self.gene_set_helper = ComboHelper(self, 'subset')

        self.gene_att_helper = ComponentIDComboHelper(self, 'gene_att',
                                                      categorical=True)

        def display_func_label(subset_group):
            return subset_group.label
        
        self.add_callback('data', self._on_data_change)
        self._on_data_change()
        
        
        self.subset_helper.choices = data_collection.subset_groups
        try:
            self.subset_helper.selection = data_collection.subset_groups[0]
        except IndexError:
            pass
        self.subset_helper.display = display_func_label
        
        
        self.gene_set_helper.choices = ['KEGG_2016','KEGG_2019_Human','KEGG_2019_Mouse','KEGG_2021_Human']
        self.gene_set_helper.selection = ['KEGG_2019_Mouse']
        

    def _on_data_change(self, *args, **kwargs):
        self.gene_att_helper.set_multiple_data([] if self.data is None else [self.data.meta['var_data']])


class DiffGeneExpState(State):

    data = SelectionCallbackProperty()
    subset1 = SelectionCallbackProperty()
    subset2 = SelectionCallbackProperty()
    #exp_att = SelectionCallbackProperty()
    gene_att = SelectionCallbackProperty()

    def __init__(self, data_collection):

        super(DiffGeneExpState, self).__init__()

        self.data_collection = data_collection
        self.data_helper = DataCollectionComboHelper(self, 'data', data_collection)
        #self.exp_att_helper = ComponentIDComboHelper(self, 'exp_att', numeric=True)
        self.gene_att_helper = ComponentIDComboHelper(self, 'gene_att',
                                                      categorical=True,
                                                      numeric=True)

        self.subset1_helper = ComboHelper(self, 'subset1')
        self.subset2_helper = ComboHelper(self, 'subset2')

        def display_func_label(subset_group):
            return subset_group.label

        self.add_callback('data', self._on_data_change)
        self._on_data_change()

        self.subset1_helper.choices = data_collection.subset_groups
        self.subset2_helper.choices = data_collection.subset_groups

        try:
            self.subset1_helper.selection = data_collection.subset_groups[0]
            self.subset2_helper.selection = data_collection.subset_groups[1]
        except IndexError:
            pass

        self.subset1_helper.display = display_func_label
        self.subset2_helper.display = display_func_label

    def _on_data_change(self, *args, **kwargs):
        #self.exp_att_helper.set_multiple_data([] if self.data is None else [self.data])
        self.gene_att_helper.set_multiple_data([] if self.data is None else [self.data.meta['var_data']])

class PCASubsetState(State):
    
    data = SelectionCallbackProperty()
    genesubset = SelectionCallbackProperty()
    do_means = CallbackProperty(False)
    do_pca = CallbackProperty(False)
    do_module = CallbackProperty(True)
    
    def __init__(self, data_collection):
        super(PCASubsetState, self).__init__()
        self.data_collection = data_collection
        self.data_helper = ManualDataComboHelper(self, 'data', data_collection)  # This should only allow cell-like datasets
        self.genesubset_helper = ComboHelper(self, 'genesubset')  # This should only allow subsets that are defined over genes...

        for data in self.data_collection:
            if data.meta.get('anndatatype') == 'obs Array':
                self.data_helper.append_data(data)


        def display_func_label(subset_group):
            return subset_group.label

        self.add_callback('data', self._on_data_change)
        self._on_data_change()
        
        self.genesubset_helper.choices = data_collection.subset_groups
        
        try:
            self.genesubset_helper.selection = data_collection.subset_groups[0]
        except IndexError:
            pass

        self.genesubset_helper.display = display_func_label
        
    def _on_data_change(self, *args, **kwargs):
        pass
        #self.gene_att_helper.set_multiple_data([] if self.data is None else [self.data.meta['var_data']])
