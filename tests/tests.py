import pytest

import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
from glue.core.component_id import ComponentID
from numpy.random import default_rng
from numpy.random import RandomState
from numpy.testing import assert_almost_equal

import anndata
import context
from glue_single_cell.data import DataAnnData
import os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print(BASE_DIR)

@pytest.fixture
def data_sparse():
    rs = RandomState(12345)
    C = rs.rand(500,700)
    C[C<0.90] = 0
    C = csc_matrix(C)
    new_adata = anndata.AnnData(X=C,dtype='float64')
    new_adata.write(filename='tests/test_data/test_dataset_float64.h5ad')
    adata = anndata.read_h5ad(os.path.join(BASE_DIR,'tests/test_data/test_dataset_float64.h5ad'),backed='r')
    d = DataAnnData(adata,label="test_name")
    yield C,d
    os.remove('tests/test_data/test_dataset_float64.h5ad')

#def test_get_sparse_data(data_sparse):
#    C,d = data_sparse
#    #print(data_sparse._components)
#    for comp in d.components:
#        print(d.get_data(comp))
    

#def test_load(data_sparse):
#    C,d = data_sparse
#    #assert type(data_sparse['X']) == anndata._core.anndata.AnnData
#    assert len(d._components)==3
    
def test_make_histogram(data_sparse):
    from scipy.sparse import find
    C,d = data_sparse
    
    xmax,ymax = C.shape
    hist_range = ((0,xmax),(0,ymax))
    bins = 11 #Histogram logic only works if range//bins = integer!! Maybe we don't care too much.
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
    #print(histogram)
    #print(correct_histogram)
    assert_almost_equal(correct_histogram, histogram)
    assert np.allclose(correct_histogram, histogram)