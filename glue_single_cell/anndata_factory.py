from glue.config import data_factory
from glue.core import Data
from pathlib import Path

from glue.core import DataCollection

from .data import DataAnnData
import anndata
import scanpy as sc

__all__ = ['is_anndata', 'read_anndata']


def df_to_data(obj,label=None):
    result = Data(label=label)
    for c in obj.columns:
        result.add_component(obj[c], str(c))
    return result

def is_anndata(filename, **kwargs):
    return filename.endswith('.h5ad') or filename.endswith('.loom')


@data_factory('AnnData data loader', is_anndata, priority=999)
def read_anndata(file_name):
    """
    Use AnnData to read a file from disk
    
    Currently supports .loom and .h5ad files
    
    
    """
    list_of_data_objs = []
    basename = Path(file_name).stem
    adata = sc.read(file_name, sparse=True, backed='r')
    
    # Get the X array as a special glue Data object
    XData = DataAnnData(adata, label=f'{basename}_X')
    list_of_data_objs.append(XData)
    
    # The var array is all components of the same length
    # and is stored by AnnData as a Pandas DataFrame
    try:
        var = adata.var
        var_data = df_to_data(var,label=f'{basename}_vars')
        var_data.add_component(adata.var_names,"var_names")
        list_of_data_objs.append(var_data)
    except:
        pass
        
    # The obs array is all components of the same length
    try:
        obs = adata.obs
        obs_data = df_to_data(obs,label=f'{basename}_obs')
        obs_data.add_component(adata.obs_names,"var_names")
        list_of_data_objs.append(obs_data)
    except:
        pass
    
    for key in adata.obsm_keys():
        
        data_arr = adata.obsm[key]
        data = Data(**{f'{key}_{i}':k for i,k in enumerate(data_arr.T)},label=f'{basename}_{key}')
        list_of_data_objs.append(data)
    
    for key in adata.varm_keys():
            
        data_arr = adata.varm[key]
        data = Data(**{f'{key}_{i}':k for i,k in enumerate(data_arr.T)},label=f'{basename}_{key}')
        list_of_data_objs.append(data)
    
    
    return list_of_data_objs