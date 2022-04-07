import os
import numpy as np
from qtpy import QtWidgets
from echo.qt import autoconnect_callbacks_to_qt

from glue.core.subset import MultiOrState
from glue.utils.qt import load_ui
from glue.core import Data, Hub, HubListener
from glue.core.message import (SubsetMessage,
                               SubsetCreateMessage,
                               SubsetUpdateMessage,
                               SubsetDeleteMessage,
                               )

from ..state import PCASubsetState
from ..anndata_factory import df_to_data

import scanpy as sc

__all__ = ['PCASubsetDialog','GeneSummaryListener']


def do_calculation_over_gene_subset(adata, genesubset, calculation = 'PCA'):
    """
    """
    print("Getting mask...")
    mask = genesubset.to_mask()
    adata_sel = adata[:, mask]  # This will fail if genesubset is not actually over genes
    if calculation == 'PCA':
        adata_sel = adata_sel.to_memory()
        sc.pp.pca(adata_sel, n_comps=10)
        data_arr = adata_sel.obsm['X_pca']
    elif calculation == 'Means':
        print("Starting mean calculation...")
        data_arr = np.expand_dims(np.sum(adata_sel.X,axis=1),axis=1)  # Expand to make same dimensionality as PCA
        print("Mean calculation finished")
    return data_arr

class GeneSummaryListener(HubListener):
    """
    A Listener to keep the new components in target_dataset
    up-to-date with any changes in the genesubset. 
    
    SubsetMessage define `subset` and `attribute` (for update?) 
    """
    def __init__(self, hub, target_dataset, genesubset, genesubset_attributes, basename, key, adata):
        self.target_dataset = target_dataset
        self.genesubset = genesubset
        self.genesubset_attributes = genesubset_attributes  #We want this to remain fixed
        self.basename = basename
        self.key = key
        self.adata = adata
        #hub.subscribe(self, SubsetCreateMessage,
        #              handler=self.update_subset)
        hub.subscribe(self, SubsetUpdateMessage,
                      handler=self.update_subset)
        hub.subscribe(self, SubsetDeleteMessage,
                      handler=self.delete_subset)
    
    def update_subset(self, message):
        """
        if the subset is the one we care about
        and if the subset still is defined over the correct attributes
        then we rerun the calculation (split that out of _apply())
        """
        subset = message.subset
        if subset == self.genesubset:
            if subset.attributes == self.genesubset_attributes:
                new_data = do_calculation_over_gene_subset(self.adata, self.genesubset, calculation = self.key)
                mapping = {f'{self.basename}_{self.key}_{i}':k for i,k in enumerate(new_data.T)}
                for x in self.target_dataset.components:  # This is to get the right component ids
                    xstr = f'{x.label}'
                    #print(xstr)
                    if xstr in mapping.keys():
                        mapping[x] = mapping.pop(xstr)
                        #del mapping[x.label]
                #print(mapping)
                #print([type(k) for k in mapping.keys()])
                self.target_dataset.update_components(mapping)
    
    def delete_subset(self, message):
        """
        Remove the attributes from target_dataset
        """
        pass
        
    def receive_message(self, message):
        print("Message received:")
        print("{0}".format(message))

class PCASubsetDialog(QtWidgets.QDialog):

    def __init__(self, collect, default=None, parent=None):

        super(PCASubsetDialog, self).__init__(parent=parent)

        self.state = PCASubsetState(collect)

        self.ui = load_ui('pca_subset.ui', self,
                          directory=os.path.dirname(__file__))
        self._connections = autoconnect_callbacks_to_qt(self.state, self.ui)

        self._collect = collect

        if default is not None:
            self.state.data = default

        self.ui.button_ok.clicked.connect(self.accept)
        self.ui.button_cancel.clicked.connect(self.reject)

    def _apply(self):
        """
        1. Take subset of anndata object (load into memory?)
        2. Run scanpy on this subset
        3. Store these values into the data somehow.
            a. If we put them into the existing AnnData object we will overwrite any existing PCAs
            b. We could return a new data object linked to the other data
            c. But the desire is to color-code e.g. a UMAP figure by these new values, which requires appending new attributes
               to the target dataset.
        
        In order to be able to make changes to the subset on-the-fly we need two things:
        1) This plug-in should establish a listener for the a specific subset and attribute
        2) We might need to make the calculation faster

        """
        genesubset = None
        
        target_dataset = self.state.data
        
        for data in self._collect:
            if target_dataset.meta['Xdata'] == data.uuid:
                Xdata = data
        
        for subset in self.state.genesubset.subsets:
            if subset.data == Xdata.meta['var_data']:  #  Find the subset on the genes, assuming we are adding to cell data
                genesubset = subset
                genesubset_attributes = subset.attributes
        if not genesubset:
            print(f"Selected subset {self.state.genesubset.label} does not seem to define genes in for {self.state.data.label}")

        adata = Xdata.Xdata
        basename = genesubset.label
        
        key = 'Means'
        data_arr = do_calculation_over_gene_subset(adata, genesubset, calculation = key)
        data = Data(**{f'{key}_{i}':k for i,k in enumerate(data_arr.T)},label=f'{basename}_{key}')
        #data.join_on_key(target_dataset,'Pixel Axis 0 [x]','Pixel Axis 0 [x]')
        #self._collect.append(data)
        for x in data.components:
            target_dataset.add_component(data.get_component(f'{x}'),f'{basename}_{x}')
        target_dataset.gene_summary_listener = GeneSummaryListener(self._collect.hub, target_dataset, genesubset, genesubset_attributes, basename, key, adata)


    @classmethod
    def summarize(cls, collect, default=None, parent=None):
        self = cls(collect, parent=parent, default=default)
        value = self.exec_()

        if value == QtWidgets.QDialog.Accepted:
            self._apply()
