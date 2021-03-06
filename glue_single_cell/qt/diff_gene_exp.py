import os
import numpy as np
from qtpy import QtWidgets
from echo.qt import autoconnect_callbacks_to_qt

from glue.core.subset import MultiOrState
from glue.utils.qt import load_ui

from ..state import DiffGeneExpState

import scanpy as sc

__all__ = ['DiffGeneExpDialog']


class DiffGeneExpDialog(QtWidgets.QDialog):

    def __init__(self, collect, default=None, parent=None):

        super(DiffGeneExpDialog, self).__init__(parent=parent)

        self.state = DiffGeneExpState(collect)

        self.ui = load_ui('diff_gene_exp.ui', self,
                          directory=os.path.dirname(__file__))
        self._connections = autoconnect_callbacks_to_qt(self.state, self.ui)

        self._collect = collect

        if default is not None:
            self.state.data = default

        self.ui.button_ok.clicked.connect(self.accept)
        self.ui.button_cancel.clicked.connect(self.reject)

    def _apply(self):
        """
        Calculate differential gene expression between the two selected
        subsets and create a new subset_group with the genes that are
        differentially expressed. 
        
        This assumes that these subsets have been defined properly on the 
        obs array
        
        Note that this copies the relevant subsets into memory as otherwise
        this won't work on an anndata gile in disk-backed mode:
        
        https://github.com/theislab/scanpy/issues/2147
        
        (the above is technically for a different scanpy function, but the same problem occurs for rank_genes_groups)
        """
        for subset in self.state.subset1.subsets:
            if subset.data == self.state.data:
                subset1 = subset
        for subset in self.state.subset2.subsets:
            if subset.data == self.state.data:
                subset2 = subset
        
        adata = self.state.data.Xdata
        
        vardata = self.state.data.meta['var_data']
        
        #mask1 = adata[subset1.label]
        #mask2 = adata[subset2.label]
        
        conditions = [
            (adata.obs[subset1.label] == '1'),
            (adata.obs[subset2.label] == '1')]
        
        choices = ['1','2']
        
        adata.obs['glue_subsets'] = np.select(conditions, choices, default='0')
         #This could be not well-defined -- but we need both masks in one var -- code mask1 as '1' and mask2 = '2' and everything else as '0'
        
        adata_selected = adata[adata.obs['glue_subsets'] != 0, :]
        try:
            adata_selected = adata_selected.to_memory()  # We should check that this is not going to be too large
        except ValueError:
            pass

        sc.tl.rank_genes_groups(adata_selected, 'glue_subsets', groups=['1'], reference='2', method='wilcoxon')
        
        n_genes = 50 #This should be a user-specified input
        gene_list = [x[0] for x in adata_selected.uns['rank_genes_groups']['names']][0:n_genes]

        print(gene_list)
        state_list = []
        for gene in gene_list:
            state_list.append(vardata.id[self.state.gene_att] == gene)
        final_state = MultiOrState(state_list)
        self.state.data_collection.new_subset_group(f'DGE between {subset1.label} and {subset2.label}', final_state)

    @classmethod
    def create_subset(cls, collect, default=None, parent=None):
        self = cls(collect, parent=parent, default=default)
        value = self.exec_()

        if value == QtWidgets.QDialog.Accepted:
            self._apply()
