from glue.config import data_factory
from glue.config import startup_action

from glue.core import Data
from glue.core import DataCollection

from glue.core.message import DataCollectionAddMessage
from glue.core import Hub, HubListener
from glue.core.qt.dialogs import warn
from glue.utils.qt import set_cursor_cm
from glue.core.qt.dialogs import warn

from qtpy import QtCore, QtWidgets
from qtpy.QtCore import Qt

from pathlib import Path
import scanpy as sc

from .data import DataAnnData

from .qt.load_data import LoadDataDialog

__all__ = ['df_to_data', 'is_anndata', 'join_anndata_on_keys', 'read_anndata', 'DataAnnDataListener', 'setup_anndata']

class AnnDataListener(HubListener):
    """
    Listen for DataAnnData objects to be added to the 
    data collection object, and, if one is, setup the
    correct join_on_key joins in a way that they will
    show up in the GUI.
    """
    def __init__(self, hub):
        hub.subscribe(self, DataCollectionAddMessage,
                      handler=self.setup_anndata)
    
    def setup_anndata(self, message):
        data = message.data
        dc = message.sender
        if isinstance(data, DataAnnData):
            setup_gui_joins(dc, data)


@startup_action("setup_anndata")
def setup_anndata(session, data_collection):
    data_collection.anndatalistener = AnnDataListener(data_collection.hub)
    return

def df_to_data(obj, label=None, skip_components=[]):
    result = Data(label=label)
    for c in obj.columns:
        if c not in skip_components:
            result.add_component(obj[c], str(c))
    return result

def is_anndata(filename, **kwargs):
    return filename.endswith('.h5ad') or filename.endswith('.loom')


def setup_gui_joins(dc, data):
    """
    Set up Join_Links that mirror the existing join_on_key links.
    This allows the user to see and remove links in the GUI.

    We cannot do this at data load because these links are defined
    at the level of a data_collection, which does not exist at
    data load time. Instead we call this through a listener
    when a DataAnnData object is added to a data collection.
    """
    try:  # If we are using a version of glue that supports links in the GUI
        from glue.core.link_helpers import JoinLink
        do_gui_link = True
    except ImportError:
        print("Cannot set up GUI join_on_key links")
        do_gui_link = False
    if do_gui_link:
        for other,joins in data._key_joins.items():
            cid, cid_other = joins
            gui_link = JoinLink(cids1=[cid[0]], cids2=[cid_other[0]], data1=data, data2=other)
            if gui_link not in dc._link_manager._external_links:
                dc.add_link(gui_link)


def join_anndata_on_keys(datasets):
    """
    Use join_on_key to stitch the various components on an anndata
    dataset back together. We join on the pixel ids and indices
    because it is far faster to do this than to join on the string
    names for genes and cells.

    TODO: Test the assumption that matching on pixel ids is really
    faster than matching on string names. Components names for 
    cells (obs) and genes (vars) are standard? in AnnData objects
    so we could join on these, which is more intuitive in the UI.
    """
    varset = {d for d in datasets if d.meta['join_on_var'] == True} 
    obsset = {d for d in datasets if d.meta['join_on_obs'] == True} 

    for dataset in datasets:
        if dataset.meta['anndatatype'] == 'X Array':
            for d in obsset:
                if (d.meta['anndatatype'] != 'X Array'):
                    dataset.join_on_key(d,'Pixel Axis 0 [y]','Pixel Axis 0 [x]')
            for d in varset:
                if (d.meta['anndatatype'] != 'X Array'):
                    dataset.join_on_key(d,'Pixel Axis 1 [x]','Pixel Axis 0 [x]')
        elif dataset.meta['anndatatype'] == 'obs Array':
            for d in obsset:
                #Do not join to self or X Array again
                if (d.meta['anndatatype'] != 'X Array') and (d.meta['anndatatype'] != 'obs Array'):
                    dataset.join_on_key(d,'Pixel Axis 0 [x]','Pixel Axis 0 [x]')
        elif dataset.meta['anndatatype'] == 'var Array':
            for d in varset:
                #Do not join to self or X Array again
                if (d.meta['anndatatype'] != 'X Array') and (d.meta['anndatatype'] != 'var Array'):
                    dataset.join_on_key(d,'Pixel Axis 0 [x]','Pixel Axis 0 [x]')
    return datasets


@data_factory("AnnData Loader", is_anndata, priority=999)
def read_anndata(file_name, skip_dialog=False, skip_components=[], subsample=False, subsample_factor=1, try_backed=False):
    """
    Use Scanpy/AnnData to read a file from disk
    
    Currently supports .loom and .h5ad files, but .loom files
    are read into memory (anndata library does not support
    a file-backed mode for them) which may cause memory issues.
    """
    list_of_data_objs = []
    basename = Path(file_name).stem

    if not skip_dialog:
        with set_cursor_cm(Qt.ArrowCursor):
            load_dialog = LoadDataDialog(filename = file_name)
            if load_dialog.exec_():
                skip_components = load_dialog.skip_components
                subsample = load_dialog.subsample
                try_backed = load_dialog.try_backed
                subsample_factor = load_dialog.subsample_factor
            else:
                return []

    if try_backed:
        try:
            adata = sc.read(file_name, sparse=True, backed='r')
            backed = True
        except OSError:
            adata = sc.read(file_name, sparse=True, backed=False)
            backed = False
    else:
        adata = sc.read(file_name, sparse=True, backed=False)
        backed = False

    if subsample:
        adata = sc.pp.subsample(adata, fraction=subsample_factor, copy=True, random_state=0)


    if backed:
        XData = DataAnnData(Xarray=adata.X, full_anndata_obj=adata, backed=backed, label=f'{basename}_X')
    else:
        XData = DataAnnData(Xarray=adata.X, backed=backed, label=f'{basename}_X')

    XData.meta['orig_filename'] = basename
    XData.meta['full_filename'] = file_name
    XData.meta['Xdata'] = XData.uuid
    XData.meta['anndatatype'] = 'X Array'
    XData.meta['join_on_obs'] = True
    XData.meta['join_on_var'] = True

    # This meta-data is attached to the DataAnnData object so that
    # We can pass it to LoadLog on save/restore
    XData.meta['loadlog_skip_components'] = skip_components
    XData.meta['loadlog_subsample'] = subsample
    XData.meta['loadlog_subsample_factor'] = subsample_factor
    XData.meta['loadlog_try_backed'] = try_backed    

    list_of_data_objs.append(XData)

    # The var array is all components of the same length
    # and is stored by AnnData as a Pandas DataFrame
    try:
        var = adata.var
        var_data = df_to_data(var, label=f'{basename}_vars', skip_components=skip_components)
        var_data.add_component(adata.var_names, "var_names")
        var_data.meta['Xdata'] = XData.uuid
        var_data.meta['anndatatype'] = 'var Array'
        var_data.meta['join_on_obs'] = False
        var_data.meta['join_on_var'] = True
        XData.meta['var_data'] = var_data 
        list_of_data_objs.append(var_data)
    except:
        pass

    for key in adata.varm_keys():
        if key not in skip_components:
            data_arr = adata.varm[key]
            data_to_add = {f'{key}_{i}':k for i,k in enumerate(data_arr.T)}
            for comp_name, comp in data_to_add.items():
                var_data.add_component(comp,comp_name)

    # The obs array is all components of the same length
    # and is stored by AnnData as a Pandas DataFrame
    try:
        obs = adata.obs
        obs_data = df_to_data(obs, label=f'{basename}_obs', skip_components=skip_components)
        obs_data.add_component(adata.obs_names, "obs_names")
        obs_data.meta['Xdata'] = XData.uuid
        obs_data.meta['anndatatype'] = 'obs Array'
        obs_data.meta['join_on_obs'] = True
        obs_data.meta['join_on_var'] = False
        XData.meta['obs_data'] = obs_data
        list_of_data_objs.append(obs_data)
    except:
        pass
    
    for key in adata.obsm_keys():
        if key not in skip_components:
            data_arr = adata.obsm[key]
            data_to_add = {f'{key}_{i}':k for i,k in enumerate(data_arr.T)}
            for comp_name, comp in data_to_add.items():
                obs_data.add_component(comp,comp_name)
    
    #obs_data.meta['xarray_data'] = Xdata
    #var_data.meta['xarray_data'] = Xdata

    return join_anndata_on_keys(list_of_data_objs)
