import os
from qtpy import QtWidgets
from echo.qt import autoconnect_callbacks_to_qt
import pandas as pd

from glue.utils.qt import load_ui

from ..state import GSEApyState

import gseapy

__all__ = ['GSEApyDialog']

def convert_genes_to_list(row):
    arr = row.split(';')
    l = [x for x in arr]
    return l

def convert_to_mouse_ids(row):
    """
    Mouse Entrez IDs are like Sox17, while enrichr returns SOX17
    This function will restore these gene IDs.
    
    Ideally we should figure out when we need to capitalize the
    gene list and only restore if we need to. 
    """
    return [x.capitalize() for x in row]

class GSEApyDialog(QtWidgets.QDialog):

    def __init__(self, collect, default=None, parent=None):

        super(GSEApyDialog, self).__init__(parent=parent)

        self.state = GSEApyState(collect)

        self.ui = load_ui('gsea.ui', self,
                          directory=os.path.dirname(__file__))
        self._connections = autoconnect_callbacks_to_qt(self.state, self.ui)

        self._collect = collect

        if default is not None:
            self.state.data = default

        self.ui.button_ok.clicked.connect(self.accept)
        self.ui.button_cancel.clicked.connect(self.reject)

    def _apply(self):
        """
        """
        if self.state.subset is not None:
            for subset in self.state.subset.subsets:
                if subset.data == self.state.data:
                    gene_subset = subset
        
        gene_list = self.state.data[self.state.gene_att] # This ignores subset for now and the need to uppercase them all
        gene_list_upper = [x for x in gene_list] #[x.upper() for x in gene_list]
        output = gseapy.enrichr(gene_list=gene_list_upper, description='pathway',
                             gene_sets=self.state.gene_set, no_plot=True)
        
        output.results['Genes'] = output.results['Genes'].apply(convert_genes_to_list)
        output.results['Genes'] = output.results['Genes'].apply(convert_to_mouse_ids)
        full_list = output.results.explode('Genes') # This is a DataFrame with each pathway-gene pair (+strength of association)
        full_list.set_index('Genes',inplace=True)
        full_term = full_list['Term'].str.get_dummies() # Expand pathways back to columns
        gene_pathways = full_term.groupby('Genes').sum() # Sum over genes associated with each expression. This is now binary
        gene_pathways.reset_index(inplace=True)
        self._collect['gene_pathways'] =  gene_pathways
        self._collect['gene_pathways'].join_on_key(self.state.data,'Genes',self.state.gene_att)
        
    @classmethod
    def enrich(cls, collect, default=None, parent=None):
        self = cls(collect, parent=parent, default=default)
        value = self.exec_()

        if value == QtWidgets.QDialog.Accepted:
            self._apply()
