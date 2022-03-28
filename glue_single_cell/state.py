from glue.core.state_objects import State
from echo import SelectionCallbackProperty
from glue.core.data_combo_helper import DataCollectionComboHelper, ComponentIDComboHelper, ComboHelper

__all__ = ['DiffGeneExpState']


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
