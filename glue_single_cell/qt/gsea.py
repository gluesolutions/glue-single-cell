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
        gene_list = gene_subset[self.state.gene_att]
        gene_list_upper = [x for x in gene_list] #[x.upper() for x in gene_list] Maybe just for humans?
        output = gseapy.enrichr(gene_list=gene_list_upper, description='pathway', organism='mouse', # Should be an option as well
                             gene_sets=self.state.gene_set, no_plot=True)
        
        self._collect[f'{self.state.gene_set} for {self.state.subset.label}'] = output.results
        
        # Theoretically we could then join on the Genes column somehow...
        # But the Genes column is Gene1;Gene2; etc. so it is not trivial and probably it is not interesting to do this.
        
        
    @classmethod
    def enrich(cls, collect, default=None, parent=None):
        self = cls(collect, parent=parent, default=default)
        value = self.exec_()

        if value == QtWidgets.QDialog.Accepted:
            self._apply()
