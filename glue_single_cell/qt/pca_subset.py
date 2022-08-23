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
from scipy.sparse import issparse
import time

__all__ = ['PCASubsetDialog','GeneSummaryListener']


def do_calculation_over_gene_subset(adata, genesubset, calculation = 'PCA'):
    """
    """
    raw = True
    #print("Getting mask...")
    mask = genesubset.to_mask()
    #Slicing the adata object is probably not the fastest thing we can do
    if calculation == 'PCA':
        adata_sel = adata[:, mask]  # This will fail if genesubset is not actually over genes
        try:
            adata_sel = adata_sel.to_memory()
        except ValueError:
            pass
        sc.pp.pca(adata_sel, n_comps=10)
        data_arr = adata_sel.obsm['X_pca']
    elif calculation == 'Module':
        adata_sel = adata[:, mask]  # This will fail if genesubset is not actually over genes

        before_memory = time.perf_counter()
        try:
            adata_sel = adata_sel.to_memory()
        except ValueError:
            pass
        after_memory = time.perf_counter()
        print(f"Loaded adata into memory in {after_memory - before_memory:0.2f} seconds", file=open('timing.txt','a'))
        gene_list = list(adata_sel.var_names)#[mask] #FIX ME We had to reapply mask IFF this was loading into memory? That seems odd
        before_score = time.perf_counter()
        try:
            sc.tl.score_genes(adata, gene_list = gene_list)
            data_arr = np.expand_dims(adata.obs['score'],axis=1)
        except ValueError:
            print("No genes found!")
            return None
        after_score = time.perf_counter()
        print(f"Calculated score in {after_score - before_score:0.2f} seconds", file=open('timing.txt','a'))

    elif calculation == 'Means':
        if raw:
            adata_sel = adata.raw.X[: , mask]  # This will fail if genesubset is not actually over genes
            if issparse(adata.raw.X):
                data_arr = np.expand_dims(adata_sel.mean(axis=1).A1,axis=1)  # Expand to make same dimensionality as PCA
        else:
            adata_sel = adata.X[: , mask]
            data_arr = np.expand_dims(adata_sel.mean(axis=1),axis=1)  # Expand to make same dimensionality as PCA
        #print("Mean calculation finished")
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
        #self.genesubset_attributes = genesubset_attributes  #We want this to remain fixed
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
            #if subset.attributes == self.genesubset_attributes:
            new_data = do_calculation_over_gene_subset(self.adata, self.genesubset, calculation = self.key)
            if new_data is not None:
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
        pass
        #print("Message received:")
        #print("{0}".format(message))

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
        1. Take a gene subset of anndata object (load into memory?)
        2. Calculate some summary statistic on this subset of genes
        3. Store these values into the data somehow.
            a. If we put them into the existing AnnData object we will overwrite any existing PCAs
            b. We could return a new data object linked to the other data
            c. But the desire is to color-code e.g. a UMAP figure by these new values, which requires appending new attributes
               to the target dataset, so this is what we do.
        
        In order to be able to make changes to the subset on-the-fly we need two things:
        1) This plug-in establishs a listener for a specific subset and attribute
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
        
        if self.state.do_means:
            key = 'Means'
        elif self.state.do_pca:
            key = 'PCA'
        elif self.state.do_module:
            key = "Module"
        data_arr = do_calculation_over_gene_subset(adata, genesubset, calculation = key)
        print(f"{data_arr.shape=}")
        if data_arr is not None:
            data = Data(**{f'{key}_{i}':k for i,k in enumerate(data_arr.T)},label=f'{basename}_{key}')
            for x in data.components:
                target_dataset.add_component(data.get_component(f'{x}'),f'{basename}_{x}')
            target_dataset.gene_summary_listener = GeneSummaryListener(self._collect.hub, target_dataset, genesubset, genesubset_attributes, basename, key, adata)


    @classmethod
    def summarize(cls, collect, default=None, parent=None):
        self = cls(collect, parent=parent, default=default)
        value = self.exec_()

        if value == QtWidgets.QDialog.Accepted:
            self._apply()
