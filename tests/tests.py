import pytest

import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
from glue.core.component_id import ComponentID
from numpy.random import default_rng
from numpy.random import RandomState
from numpy.testing import assert_almost_equal, assert_equal
from scipy.sparse import find

import anndata
import context
from glue_single_cell.data import DataAnnData
import os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print(BASE_DIR)

SPARSE_BACKED_OBS_NUM = 500
SPARSE_BACKED_VAR_NUM = 700

def print_size_in_MB(x):
    """
    Convenience method for computing the size of objects
    """
    print('{:.3} MB'.format(x.__sizeof__()/1e6))

@pytest.fixture
def data_sparse_backed():
    rs = RandomState(12345)
    C = rs.rand(SPARSE_BACKED_OBS_NUM,SPARSE_BACKED_VAR_NUM)
    C[C<0.90] = 0
    C = csc_matrix(C)
    new_adata = anndata.AnnData(X=C,dtype='float64')
    new_adata.write(filename='tests/test_data/test_dataset_float64.h5ad')
    adata = anndata.read_h5ad(os.path.join(BASE_DIR,'tests/test_data/test_dataset_float64.h5ad'),backed='r')
    d = DataAnnData(adata,label="test_name")
    yield C,d
    os.remove('tests/test_data/test_dataset_float64.h5ad')

def test_get_data_view_sparse_backed(data_sparse_backed):
    C,d = data_sparse_backed
    views = (np.s_[:],
             np.s_[0,:],
             np.s_[:,0],
             np.s_[:,10:30],
             np.s_[37:87,:],
             )
    for view in views:
        for comp in d.components:
            if comp in d._pixel_component_ids:
                pass
            else:
                comp_data = d.get_data(comp,view)
                #print(comp_data)
                #print(comp_data.shape)
                #print(C[view].shape)
                assert comp_data.size == C[view].size
                assert_equal(comp_data,C[view].data)


def test_get_data_sparse_backed(data_sparse_backed):
    C,d = data_sparse_backed
    #print(data_sparse._components)
    for comp in d.components:
        comp_data = d.get_data(comp)
        if comp in d._pixel_component_ids:
            pass
        else:
            assert len(comp_data) == C.size



def test_data_setup_sparse_backed(data_sparse_backed):
    C,d = data_sparse_backed
    #assert type(data_sparse['X']) == anndata._core.anndata.AnnData
    assert len(d._components)==3
    
def test_make_histogram_sparse_backed(data_sparse_backed):
    C,d = data_sparse_backed
    xmax,ymax = C.shape
    hist_range = ((0,xmax),(0,ymax))
    bins = (11,11) #Histogram logic only works if range//bins = integer!! Maybe we don't care too much.
    x,y,w = find(C)
    correct_histogram = np.histogram2d(x = x,
                                       y = y,
                                       bins=bins,
                                       range=hist_range,
                                       weights = w)[0]
    histogram = d.compute_histogram(cids=['Pixel Axis 0 [y]','Pixel Axis 1 [x]'],
                                              bins=bins,
                                              range=hist_range,
                                              weights = ['X'])
    assert_almost_equal(correct_histogram, histogram)
    assert np.allclose(correct_histogram, histogram)