import os
import pytest
import numpy as np
from numpy.testing import assert_array_equal
from glue.core.data_collection import DataCollection
from glue.core import Data
from glue.app.qt import GlueApplication
from unittest.mock import patch
from glue.core.state import GlueUnSerializer

from ..gsea import GSEApyDialog
from ..pca_subset import dialog, PCASubsetDialog

class TestEnrichr(object):


    def get_data(self, try_backed=False):

        fake_data = Data(gene_id = ['Sox17','Adhfe1','Prex2','Msc','Rdh10','Gsta3'],
                   qtl = [1,2,3,4,5,5],
                   label = 'qtl')
        return fake_data

    def setup_method(self):
        self.app = GlueApplication()
        self.session = self.app.session
        self.hub = self.session.hub
    
        self.dc = self.session.data_collection
        self.dc.append(self.get_data())

        self.gene_data = self.dc[0]        

        self.subset1 = self.dc.new_subset_group(label='Interesting Genes', subset_state=self.gene_data.id['qtl'] < 5)

    def do_enrichr(self):
        gseadiag = GSEApyDialog(self.dc)
        gseadiag.state.data = self.dc[0]
        gseadiag.state.subset = self.subset1
        gseadiag.state.gene_att = self.dc[0].id['gene_id']
        gseadiag.state.gene_set = 'KEGG_2019_Mouse'
    
        gseadiag._apply(do_dialog=False)


    def test_get_enrich(self):
        assert len(self.dc) == 1 # We add a new dataset
        self.do_enrichr()
        assert len(self.dc) == 2 # We add a new dataset
        assert len(self.dc[1].components) == 11

    def test_save_and_restore(self, tmpdir):
        self.do_enrichr()
        filename = tmpdir.join('test_anndata_load_session.glu').strpath
        self.session.application.save_session(filename)
        with open(filename, 'r') as f:
            session = f.read()
    
        state = GlueUnSerializer.loads(session)
    
        ga = state.object('__main__')
    
        dc = ga.session.data_collection
        assert len(dc) == 2 # Check that the new dataset save/restores