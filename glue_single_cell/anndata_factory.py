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


def join_anndata_on_keys(datasets):
    """
    Use join_on_key to stitch the various components on an anndata
    dataset back together. We join on the pixel ids and indices
    because it is far faster to do this than to join on the string
    names for genes and cells.
    """
    varset = {d for d in datasets if d.meta['join_on_var'] == True} 
    obsset = {d for d in datasets if d.meta['join_on_obs'] == True} 
    
    for dataset in datasets:
        if dataset.meta['anndatatype'] == 'X Array':
            for d in varset:
                if d.meta['anndatatype'] != 'X Array': #Do not join to self
                    dataset.join_on_key(d,'Pixel Axis 0 [y]','Pixel Axis 0 [x]')
            for d in obsset:
                if d.meta['anndatatype'] != 'X Array': #Do not join to self
                    dataset.join_on_key(d,'Pixel Axis 1 [x]','Pixel Axis 0 [x]')
        elif dataset.meta['anndatatype'] == 'obsm Array':
            for d in obsset:
                #Do not join to self or X Array again
                if (d.meta['anndatatype'] != 'X Array') and (d.meta['anndatatype'] != 'obsm Array'):
                    dataset.join_on_key(d,'Pixel Axis 0 [x]','Pixel Axis 0 [x]')
        elif dataset.meta['anndatatype'] == 'varm Array':
            for d in varset:
                #Do not join to self or X Array again
                if (d.meta['anndatatype'] != 'X Array') and (d.meta['anndatatype'] != 'varm Array'):
                    dataset.join_on_key(d,'Pixel Axis 0 [x]','Pixel Axis 0 [x]')
    return datasets



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
    XData.meta['anndatatype'] = 'X Array'
    XData.meta['join_on_obs'] = True
    XData.meta['join_on_var'] = True
    
    
    list_of_data_objs.append(XData)
    
    # The var array is all components of the same length
    # and is stored by AnnData as a Pandas DataFrame
    try:
        var = adata.var
        var_data = df_to_data(var,label=f'{basename}_vars')
        var_data.add_component(adata.var_names,"var_names")
        var_data.meta['anndatatype'] = 'var Array'
        var_data.meta['join_on_obs'] = False
        var_data.meta['join_on_var'] = True
    
    
        list_of_data_objs.append(var_data)
    except:
        pass
        
    # The obs array is all components of the same length
    try:
        obs = adata.obs
        obs_data = df_to_data(obs,label=f'{basename}_obs')
        obs_data.add_component(adata.obs_names,"var_names")
        obs_data.meta['anndatatype'] = 'obs Array'
        obs_data.meta['join_on_obs'] = True
        obs_data.meta['join_on_var'] = False
    
    
        list_of_data_objs.append(obs_data)
    except:
        pass
    
    for key in adata.obsm_keys():
        
        data_arr = adata.obsm[key]
        data = Data(**{f'{key}_{i}':k for i,k in enumerate(data_arr.T)},label=f'{basename}_{key}')
        data.meta['anndatatype'] = 'obsm Array'
        data.meta['join_on_obs'] = True
        data.meta['join_on_var'] = False
    
    
        list_of_data_objs.append(data)
    
    for key in adata.varm_keys():
            
        data_arr = adata.varm[key]
        data = Data(**{f'{key}_{i}':k for i,k in enumerate(data_arr.T)},label=f'{basename}_{key}')
        data.meta['anndatatype'] = 'varm Array'
        data.meta['join_on_var'] = True
        data.meta['join_on_obs'] = False
        list_of_data_objs.append(data)
    
    join_anndata_on_keys(list_of_data_objs)
    return list_of_data_objs