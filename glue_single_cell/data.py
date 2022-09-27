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

from collections import OrderedDict
import uuid

import pandas as pd
import numpy as np

from glue.core.component import Component, CoordinateComponent, DerivedComponent
from glue.core.data import BaseCartesianData, Data
from glue.core.component_id import ComponentID, ComponentIDDict
from glue.core.component_id import ComponentIDList

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
    definitions
    
    SubsetMessage define `subset` and `attribute` (for update?) 
    """
    def __init__(self, hub, anndata):
        self.anndata = anndata
        self.Xdata = anndata.Xdata
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
                self.Xdata.obs[subset.label] = obs_mask.astype('int').astype('str')
            except IncompatibleAttribute:
                pass
        elif subset.data == self.anndata.meta['var_data']:
            try:
                var_mask = subset.to_mask()
                self.Xdata.var[subset.label] = var_mask.astype('int').astype('str')
            except IncompatibleAttribute:
                pass

    def delete_subset(self, message):
        if subset.data == self.anndata.meta['obs_data']:
            try:
                self.Xdata.obs.drop(subset.label,axis=1,inplace=True)
            except:
                pass
        elif subset.data == self.anndata.meta['var_data']:
            try:
                self.Xdata.var.drop(subset.label,axis=1,inplace=True)
            except:
                pass


    def receive_message(self, message):
        pass
        #print("Message received:")
        #print("{0}".format(message))

class AnnData(Data):
    def __init__(self, data, label='', **kwargs):
        super().__init__(label=label, **kwargs)
        component_id = ComponentID(label='X', parent=self)
        comp_to_add = Component(data.X)
        self.add_component(comp_to_add,label=component_id)
        self._components[component_id] = comp_to_add
        self._shape = data.shape
        print(f"{self._shape=}")
        print(f"{comp_to_add.shape=}")

    def attach_subset_listener(self):
        if self.hub is not None:
            self.subset_listener = SubsetListener(self.hub, self)


class DataAnnData(Data):
    
    def __init__(self, data, label="", coords=None, **kwargs):
        """
        data is an AnnData object with a single data matrix of
        shape #observations x #variables. This can be either
        a dense numpy array or a scipy sparse array and can
        be either in memory or (more typically) on disk.
        
        The other components of an AnnData object:
        var, obs, varm, obsm are stored as separate (regular)
        glue data objects, because they have different
        dimensions. They are connected to this data object
        through join_on_key relations. See anndata_factory.py
        for how this is done in the loader.
        """
        super(BaseCartesianData, self).__init__()
        #We sanitize the underlying data array because
        #A) Currently, many scanpy functions call sanitize_anndata()
        #B) sanitize does not work on a slice/subset, so we do it here
        self.Xdata = data
        scanpy._utils.sanitize_anndata(self.Xdata)
        self.label = label

        self._shape = ()

        # Components
        self._components = OrderedDict()
        self._externally_derivable_components = OrderedDict()
        self._pixel_aligned_data = OrderedDict()
        self._pixel_component_ids = ComponentIDList()
        self._world_component_ids = ComponentIDList()

        # Coordinate conversion object
        self.coords = coords

        self.id = ComponentIDDict(self)

        self._coordinate_links = []

        self.edit_subset = None

        self._key_joins = {}

        component_id = ComponentID(label='X', parent=self)
        comp_to_add = Component(self.Xdata)
        self.add_component(comp_to_add,label=component_id)
        self._components[component_id] = comp_to_add
        self._shape = self.Xdata.shape
        self._label = label
        # To avoid circular references when saving objects with references to
        # the data, we make sure that all Data objects have a UUID that can
        # uniquely identify them.
        self.uuid = str(uuid.uuid4())
        
        # When we create this data object, we don't have a hub set up yet,
        # so we can't init the Listener at Data creation time. Instead,
        # we add a Listener to the DataCollection object in a custom
        # startup_action, and this adds a subset_listener to DataAnnData
        # objects that are added to the data collection.
                
    def attach_subset_listener(self):
        if self.hub is not None:
            self.subset_listener = SubsetListener(self.hub, self)
        
        
    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, value):
        if getattr(self, '_label', None) != value:
            self._label = value
            self.broadcast(attribute='label')
        elif value is None:
            self._label = value

    @property
    def shape(self):
        return self._shape

    @property
    def main_components(self):
        return [c for c in self.component_ids() if
                not isinstance(self._components[c], (DerivedComponent, CoordinateComponent))]
    
    def get_kind(self, cid):
    
        comp = self.get_component(cid)
        
        try:
            dtype = comp.dtype
        except AttributeError: #This will happen for AnnData component
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

    def get_mask(self, subset_state, view=None):
        return subset_state.to_mask(self,view=view) #Is this sufficient?

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

# Defining something like this would allow us to get the Xdata object from our DataAnnData object
# more transparently
# @data_translator(anndata.AnnData)
# class GeoPandasTranslator:
#  
#     def to_data(self, data):
#         return GeoRegionData(data)
#  
#     def to_object(self, data_or_subset, attribute=None):
#         gdf = geopandas.GeoDataFrame()
#         coords = data_or_subset.coordinate_components
#         if isinstance(data_or_subset, Subset):
#             #geom = data_or_subset.data.geometry
#             centroids = data_or_subset.data._centroid_component_ids #because these are sort of fake coords
#             crs = data_or_subset.data.meta['crs']
#         else:
#             #geom = data_or_subset.geometry
#             centroids = data_or_subset._centroid_component_ids #because these are sort of fake coords
#             crs = data_or_subset.meta['crs']
#             
#         #gdf.geometry = geom
#         for cid in data_or_subset.components:
#             if (cid not in coords) and (cid not in centroids):
#                 if cid.label == 'geometry':
#                     g = geopandas.GeoSeries.from_wkt(data_or_subset[cid])
#                     gdf[cid.label] = g
#                 else:
#                     gdf[cid.label] = data_or_subset[cid]
#         gdf.set_geometry("geometry",inplace=True)
#         gdf.crs = crs
#         return gdf
