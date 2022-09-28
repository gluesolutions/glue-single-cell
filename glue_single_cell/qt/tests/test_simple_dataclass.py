import numpy as np
import os
from glue.core.component import Component
from glue.core.component_id import ComponentID
import pytest
from glue.core.data import BaseCartesianData, Data
from glue.utils import view_shape
from glue.config import data_factory

from glue.app.qt import GlueApplication
from glue.core.state import (GlueSerializer, GlueUnSerializer,
                     saver, loader, VersionedDict)
import anndata as ad
from scipy.sparse import csr_matrix
import scanpy as sc
from glue_single_cell.anndata_factory import read_anndata, read_anndata_new

DATA = os.path.join(os.path.dirname(__file__), 'data')

@saver(csr_matrix)
def _save_csr(csr, context):
    return {'sparse_matrix':"sparse_matrix"}

@loader(csr_matrix)
def _load_csr(rec, context):
    pass


class AnnDataMock(Data):

    def __init__(self, **kwargs):
        super().__init__()

        #counts = csr_matrix(np.random.poisson(1, size=(100, 2000)), dtype=np.float32)
        adata = ad.read(os.path.join(DATA, 'test_data.h5ad'))

        component_id = ComponentID(label='X', parent=self)
        comp_to_add = Component(adata.X)
        self.add_component(comp_to_add,label=component_id)
        self._components[component_id] = comp_to_add
        self._shape = adata.shape
        print(f"{self._shape=}")
        print(f"{comp_to_add.shape=}")


if __name__ == "__main__":
    app = GlueApplication()
    session = app.session
    data_collection = session.data_collection
    app.load_data(os.path.join(DATA, 'test_data.h5ad'), factory=read_anndata, skip_dialog=True)
    #data = Data(x=[1,2,3],label='test')
    #data = AnnData()
    #import ipdb; ipdb.set_trace()
    #data_collection.append(data)
    filename = 'test_read_ann_data_thing.glu'
    app.save_session(filename, include_data=False)


def test_simple_create_save_restore(tmpdir):
    app = GlueApplication()
    session = app.session
    data_collection = session.data_collection
    data = AnnData()
    data_collection.append(data)
    
    #filename = 'yo_delete.glu'
    filename = tmpdir.join('test_anndata_load_session.glu').strpath
    app.save_session(filename, include_data=False)
    
    with open(filename, 'r') as f:
        session = f.read()
    app2 = GlueApplication.restore_session(filename)
    
    dc = app2.session.data_collection
    assert len(dc) == 1
    ga.close()