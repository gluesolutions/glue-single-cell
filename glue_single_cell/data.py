"""
A glue Dataclass that wraps an AnnData object

The primary motivation is to support on-disk access
to a dataset that (even sparse) may be too large
to fit comfortably in memory.

We store this data in a dual format, as both a native AnnData
object AND as glue components. This allows us to call scanpy
functions on the AnnData object without much additional
bookkeeping. 

AnnData objects include many things of different dimensions.
This DataAnnData class only exposes the X matrix of data values
with dimension num_obs x num_vars as a glue component. The other
parts of the AnnData object are stored as regular glue data
objects -- their creation and linking with the DataAnnData is
handled by the data loader.

All the obsm arrays and obs table are combined as one dataset
All the varm arrays and var table are combined as one dataset

The obsp and varp arrays could be extremely large and require
a dedicated data class. We do not deal with these yet.

Because of automatic sanitization 
(see __[this issue](https://github.com/theislab/scanpy/issues/1747)__)
we generally have to sanitize an AnnData object before using it in
Scanpy. Since sanitization cannot be done on a subset, we sanitize
the full dataset at initialization time.

TODO: The anndata/scanpy convention of updating the data
structure when new calculations are performed only works
with glue if we know how to add components to the appropriate
glue datasets (which currently are only linked-by-key) in a 
loose way, and should probably be more tightly coupled. 

anndata - Annotated Data
https://anndata.readthedocs.io/en/latest/

Scanpy - Single-Cell Analysis in Python
https://scanpy.readthedocs.io/en/stable/
"""

import uuid

import pandas as pd
import numpy as np
from scipy.sparse import isspmatrix
from collections import OrderedDict

from glue.core.component import Component, CoordinateComponent, DerivedComponent
from glue.core.data import BaseCartesianData, Data, Subset
from glue.core.component_id import ComponentID, ComponentIDDict, PixelComponentID, ComponentIDList

from glue.core.component_link import ComponentLink, CoordinateComponentLink
from glue.core.exceptions import IncompatibleAttribute
    
from fast_histogram import histogram1d, histogram2d

from glue.core import data_factories as df
from glue.core import Hub, HubListener

from glue.core.message import (DataMessage,
                                DataCollectionMessage,
                                SubsetMessage,
                                SubsetCreateMessage,
                                SubsetUpdateMessage,
                                SubsetDeleteMessage,
                               )
from glue.core.state import (GlueSerializer, GlueUnSerializer,
                                saver, loader, VersionedDict)
from glue.config import session_patch, data_translator

import anndata
import scanpy




def get_subset(subset_name, data_collection, app=None, save_to_disk=False):
    """
    Return a view of the anndata object that corresponds
    to the desired subset
    """
    slice_list,uuid = get_subset_mask(subset_name, data_collection)
    
    for dataset in data_collection:
        try:
            data_uuid = dataset.meta['Xdata']
            orig_filename = dataset.meta['orig_filename']
        except KeyError:
            return None
        if isinstance(dataset, DataAnnData) and (data_uuid == uuid):
            goodsubset = None
            try:
                goodsubset = dataset.Xdata[slice_list,:]
            except IndexError:
                goodsubset = dataset.Xdata[:,slice_list]
            if goodsubset:
                if save_to_disk:
                    new_filename = f'{orig_filename}_{subset_name.replace(" ","")}.h5ad'
                    newadata = goodsubset.copy(filename=new_filename) #This creates the file but leaves if open in the original mode
                    newadata.file.close() #So we close it
                    goodsubset = app.load_data(new_filename)[0].Xdata #This adds the data to the data_collection and returns a reference to the first/main dataset
                    #The problem with this approach is that we keep adding new data
                    #to the data collection every time we update the subset and call
                    #this function again. We could... delete from the data object and force a reload?
                    #What does this look like in the case were we do not save to disk?
                    
                    
                    #goodsubset = anndata.read(new_filename,backed='r+') #And reopen it so that it is editable
                    print(f"Subset is being written to disk as {new_filename}")                
                else:
                    print(f"Subset is defined as a slice of the full dataset at {orig_filename}. To save this subset as a new object to disk, re-run this command with save_to_disk=True")    
                return goodsubset
                
def get_subset_mask(subset_name, data_collection):
    """
    Get a subset mask for a given subset_name by
    trying to get the mask for all datasets in the
    data_collection.
    """
    target_subset = None
    for subset in data_collection.subset_groups:
        if subset.label == subset_name:
            target_subset = subset
    if target_subset == None:
        print(f"Subset {subset_name} not found. Is this subset defined?")
        return None
    for subset_on_data in target_subset.subsets:
        try:
            list_of_obs = subset_on_data.to_mask()
            try:
                uuid = subset_on_data.data.meta['Xdata']
            except KeyError:
                uuid = None
            return list_of_obs, uuid
        except IncompatibleAttribute:
            pass


class SubsetListener(HubListener):
    """
    A Listener to keep the X array of a DataAnnData
    object updated with the current glue subset
    definitions. We use this for the keeping the
    Xdata array up-to-date so that we can use it
    in some of the plug-ins.
    
    This seems like a fairly expensive option though
    We should only create an AnnData object when we 
    need it
    
    """
    def __init__(self, hub, anndata):
        self.anndata = anndata
        #self.Xdata = anndata.Xdata
        hub.subscribe(self, SubsetCreateMessage,
                      handler=self.update_subset)
        hub.subscribe(self, SubsetUpdateMessage,
                      handler=self.update_subset)
        hub.subscribe(self, SubsetDeleteMessage,
                      handler=self.delete_subset)

    def update_subset(self, message):
        """
        The trick here is that we want to 
        react to subsets on the obs and var arrays
        """
        
        subset = message.subset
        if subset.data == self.anndata.meta['obs_data']:
            try:
                obs_mask = subset.to_mask()
                #self.Xdata.obs[subset.label] = obs_mask.astype('int').astype('str')
            except IncompatibleAttribute:
                pass
        elif subset.data == self.anndata.meta['var_data']:
            try:
                var_mask = subset.to_mask()
                #self.Xdata.var[subset.label] = var_mask.astype('int').astype('str')
            except IncompatibleAttribute:
                pass

    def delete_subset(self, message):
        if subset.data == self.anndata.meta['obs_data']:
            try:
                pass
                #self.Xdata.obs.drop(subset.label,axis=1,inplace=True)
            except:
                pass
        elif subset.data == self.anndata.meta['var_data']:
            try:
                pass
                #self.Xdata.var.drop(subset.label,axis=1,inplace=True)
            except:
                pass


    def receive_message(self, message):
        pass
        #print("Message received:")
        #print("{0}".format(message))

class DataAnnData(Data):
    def __init__(self, label='', full_anndata_obj=None, backed=False, coords=None, **kwargs):
        super().__init__(label=label, coords=None)

        self._Xdata = full_anndata_obj
        self.backed = backed
        self.sparse = False

        # The normal add_component in Data.__init__ fails for a sparse array
        # because it tries to use autotyped().
        for lbl, Xdata in sorted(kwargs.items()):
            component_id = ComponentID(label=lbl, parent=self)
            comp_to_add = Component(Xdata)
            self.add_component(comp_to_add,label=component_id)
            if isspmatrix(Xdata):
                self.sparse = True
            elif type(Xdata) == anndata._core.sparse_dataset.SparseDataset:
                self.sparse = True

        # When we create this data object, we don't have a hub set up yet,
        # so we can't init the SubsetListener at Data creation time. Instead,
        # we add a Listener to the DataCollection object in a custom
        # startup_action, and this adds a subset_listener to DataAnnData
        # objects that are added to the data collection.

    @property
    def Xdata(self):
        if self._Xdata:
            return self._Xdata
        else:
            return self.get_object(cls=anndata.AnnData)

    @Xdata.setter
    def Xdata(self, value):
        self._Xdata = value

    def get_mask(self, subset_state, view=None):
        """
        This is a bit of a hack. Because we have join_on_keys on both
        dimensions of the X array, Data.get_mask() will normally
        default to trying get_mask_with_key_joins. By not 
        allowing that we prevent subsets on obs/var arrays traveling
        through the X array to completely cover var/ob arrays.
        """
        return subset_state.to_mask(self,view=view)
        
    def attach_subset_listener(self):
        if self.hub is not None:
            self.subset_listener = SubsetListener(self.hub, self)

    def get_kind(self, cid):

        comp = self.get_component(cid)

        try:
            dtype = comp.dtype
        except AttributeError: #This might happen on the X array 
            return 'numeric'

        if comp.datetime:
            return 'datetime'
        elif comp.numeric:
            return 'numerical'
        elif comp.categorical:
            return 'categorical'
        else:
            raise TypeError("Unknown data kind")

    def get_data(self, cid, view=None):
        """
        This is the tricky function for anndata backed store, since
        in general returning the full data object leads to memory problems.
        
        We also might have an issue that anndata does not support views
        into other views, so we need to decompose things.
        
        We could have get_data return an iterator directly. 
        """

        if not self.backed:
            return super().get_data(cid, view=view)

        if isinstance(cid, ComponentLink):
            return cid.compute(self, view)

        if cid in self._components:
            comp = self._components[cid]
        elif cid in self._externally_derivable_components:
            comp = self._externally_derivable_components[cid]
        else:
            raise IncompatibleAttribute(cid)

        if view is not None:
            result = comp[view]
            if cid not in self._pixel_component_ids:
                result = self.Xdata.X[view].data #This probably just loads everything into memory
        else:
            if cid not in self._pixel_component_ids:
                result = self.Xdata.X[:,:].data #This loads everything into memory
            else:
                result = comp.data

        return result

    def compute_statistic(self, statistic, cid,
                          axis=None, finite=True,
                          positive=False, subset_state=None,
                          percentile=None, random_subset=None):
        #We can probably just convert subset_state to a mask and pass it as a view?
        if statistic == 'minimum':
            #self.data[:,:].X.min(axis=axis) #Another way to do this
            return np.min(self.Xdata.X,axis=axis) #Some concern this might not be correct?
        elif statistic == 'maximum':
            return np.max(self.Xdata.X,axis=axis)
        elif statistic == 'mean':
            return np.mean(self.Xdata.X,axis=axis)
        elif statistic == 'median':
            return np.median(self.Xdata.X,axis=axis) #We can't do the other way for this
        elif statistic == 'percentile':
            return np.percentile(self.Xdata.X,percentile,axis=axis) #We can't do the other way for this
        elif statistic == 'sum':
            return np.sum(self.Xdata.X,axis=axis) #We can't do the other way for this
        
    def compute_histogram(self, cids, weights=None, range=None, bins=None, log=None, subset_state=None, chunk_size=1000):
        """
        Compute an n-dimensional histogram with regularly spaced bins.

        We restrict this to certain cids in order to enable fast histogram
        calculations.

        Parameters
        ----------
        cids : list of str or `ComponentID`
            Component IDs to compute the histogram over
        weights : str or ComponentID
            Component IDs to use for the histogram weights
        range : list of tuple
            The ``(min, max)`` of the histogram range
        bins : list of int
            The number of bins
        log : list of bool
            Whether to compute the histogram in log space
        subset_state : `SubsetState`, optional
            If specified, the histogram will only take into account values in
            the subset state.
        chunk_size : number of rows per chunk for calculating histogram
        """
        from scipy.sparse import find
        
        #TODO: We should return a NotImplementedError if we try to do a histogram
        #on anything that this not the gene/cell coords and X matrix
        
        if len(cids) > 2:
            raise NotImplementedError()
        
        ndim = len(cids)
        
        if ndim == 1:
            xmin, xmax = range[0]
            xmin, xmax = sorted((xmin, xmax))
        else:
            (xmin, xmax), (ymin, ymax) = range
            xmin, xmax = sorted((xmin, xmax))
            ymin, ymax = sorted((ymin, ymax))
        
        if ndim >= 1:
            xmax += 10 * np.spacing(xmax)
        if ndim >= 2:
            ymax += 10 * np.spacing(ymax)
                
        if ndim == 1:
            range = (xmin, xmax)
            chunked_histogram = np.zeros(bins[0])
            for chunk, start, end in self.Xdata.chunked_X(chunk_size):
                x,y,w=find(chunk)
                chunked_histogram += histogram1d(x+start, bins=bins[0], range=range, weights = w)
            return chunked_histogram
            
        if ndim > 1:
            range = [(xmin, xmax), (ymin, ymax)]
            chunked_histogram = np.zeros(bins)
            for chunk, start, end in self.Xdata.chunked_X(chunk_size):
                x,y,w=find(chunk)
                chunked_histogram += np.histogram2d(x+start, y, bins=bins, range=range, weights = w)[0]
            return chunked_histogram

@saver(DataAnnData, version=1)
def _save_anndata(data, context):
    """
    Custom save function for DataAnnData.
    Much of this duplicates _save_data_5 for Data
    We probably do not actually *need* this
    """
    result = dict(components=[(context.id(c),
                             context.id(data.get_component(c)))
                             for c in data._components],
                 subsets=[context.id(s) for s in data.subsets],
                 label=data.label)

    if data.coords is not None:
        result['coords'] = context.id(data.coords)

    result['style'] = context.do(data.style)

    def save_cid_tuple(cids):
        return tuple(context.id(cid) for cid in cids)

    result['_key_joins'] = [[context.id(k), save_cid_tuple(v0), save_cid_tuple(v1)]
                            for k, (v0, v1) in data._key_joins.items()]
    result['uuid'] = data.uuid

    result['primary_owner'] = [context.id(cid) for cid in data.components if cid.parent is data]
    # Filter out keys/values that can't be serialized
    meta_filtered = OrderedDict()
    for key, value in data.meta.items():
        try:
            context.do(key)
            context.do(value)
        except GlueSerializeError:
            continue
        else:
            meta_filtered[key] = value
    result['meta'] = context.do(meta_filtered)
    return result


@session_patch(priority=0)
def correct_coords_problem(rec):
    """
    The default LoadLog incorrectly sets force_coords = True
    for the anndata_factory (because of an unfortunate
    coincidence where there are four coordinate components in
    the full list of coordinate components and X.ndim = 2). 

    We patch the LoadLog here never force_coords
    """
    for key, value in rec.items():
        if key=='LoadLog':
            value['force_coords'] = False            


@loader(DataAnnData, version=1)
def _load_anndata(rec, context):
    """
    Custom load function for AnnData.
    This is the same as the chain of logic in 
    _save_data_5 for Data, but result is an AnnData object
    instead.
    """

    label = rec['label']
    result = DataAnnData(label=label)

    # we manually rebuild pixel/world components, so
    # we override this function. This is pretty ugly
    result._create_pixel_and_world_components = lambda ndim: None

    comps = [list(map(context.object, [cid, comp]))
             for cid, comp in rec['components']]

    for icomp, (cid, comp) in enumerate(comps):
        if isinstance(comp, CoordinateComponent):
            comp._data = result

            # For backward compatibility, we need to check for cases where
            # the component ID for the pixel components was not a PixelComponentID
            # and upgrade it to one. This can be removed once we no longer
            # support pre-v0.8 session files.
            if not comp.world and not isinstance(cid, PixelComponentID):
                cid = PixelComponentID(comp.axis, cid.label, parent=cid.parent)
                comps[icomp] = (cid, comp)

        result.add_component(comp, cid)

    assert result._world_component_ids == []

    coord = [c for c in comps if isinstance(c[1], CoordinateComponent)]
    coord = [x[0] for x in sorted(coord, key=lambda x: x[1])]

    if getattr(result, 'coords') is not None:
        assert len(coord) == result.ndim * 2
        result._world_component_ids = coord[:len(coord) // 2]
        result._pixel_component_ids = coord[len(coord) // 2:]
    else:
        assert len(coord) == result.ndim
        result._pixel_component_ids = coord

    # We can now re-generate the coordinate links
    result._set_up_coordinate_component_links(result.ndim)

    for s in rec['subsets']:
        result.add_subset(context.object(s))    

    result.style = context.object(rec['style'])

    if 'primary_owner' in rec:
        for cid in rec['primary_owner']:
            cid = context.object(cid)
            cid.parent = result
    yield result

    def load_cid_tuple(cids):
        return tuple(context.object(cid) for cid in cids)

    result._key_joins = dict((context.object(k), (load_cid_tuple(v0), load_cid_tuple(v1)))
                             for k, v0, v1 in rec['_key_joins'])
    if 'uuid' in rec and rec['uuid'] is not None:
        result.uuid = rec['uuid']
    else:
        result.uuid = str(uuid.uuid4())
    if 'meta' in rec:
        result.meta.update(context.object(rec['meta']))


@data_translator(anndata.AnnData)
class DataAnnDataTranslator:

    def to_data(self, data):
        """
        This is not quite right. We actually
        need to return a bunch of linked datasets
        (all the logic in read_anndata)
        """
        return DataAnnData(data.X) 

    def unwrap_components(self, data):
        output_dict = {}
        for cid in data.components:
            if not isinstance(cid, PixelComponentID):
                output_dict[cid.label] = data[cid]
            if isinstance(cid, PixelComponentID):
                output_index = data[cid]
        output_df = pd.DataFrame(output_dict, index=output_index)

    def to_object(self, data_or_subset, attribute=None):
        """
        Re-create an AnnData object from the Glue data objects
        This does not faithfully reproduce things since varm and
        obsm arrays have been flattened with var and obs arrays
        and some things (like .uns and layers) are ignored

        How does this work with backed mode?
        """
     
        if isinstance(data_or_subset, Subset):
            pass
        elif isinstance(data_or_subset, DataAnnData):
            
            obs_glue_data = data_or_subset.meta['obs_data']
            var_glue_data = data_or_subset.meta['var_data']
            obs_data = self.unwrap_components(obs_glue_data)
            var_data = self.unwrap_components(var_glue_data)

            adata = anndata.AnnData(data_or_subset['Xarray'], obs=obs_data, var=var_data)
            return adata
        else:
            pass
            