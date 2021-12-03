"""
A glue dataclass that wraps an AnnData object

The primary motivation is to support on-disk access
to a sparse datatype

There is a lot in Data that will probably just work normally
here. Maybe we should sub-class that instead?

"""

from collections import OrderedDict
import uuid

import pandas as pd
import numpy as np

from glue.core.component import Component
from glue.core.data import BaseCartesianData, Data
from glue.core.component_id import ComponentID, ComponentIDDict
from glue.core.component_id import ComponentIDList

from glue.utils import (compute_statistic, unbroadcast, iterate_chunks,
    datetime64_to_mpl, broadcast_to, categorical_ndarray,
    format_choices, random_views_for_dask_array)
    
from fast_histogram import histogram1d, histogram2d

import anndata

class DataAnnData(Data):
    def __init__(self, data, label="", coords=None, **kwargs):
        """
        data is an AnnData object with a single sparse
        matrix and (generally) with a filename backing
        """
        super(BaseCartesianData, self).__init__()
        #print(data)
        self.Xdata = data
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

        #for lbl, data in sorted(kwargs.items()):
        #    self.add_component(data, lbl)

        self._key_joins = {}

        #self.add_component(data,'X')
        component_id = ComponentID(label='X', parent=self)
        comp_to_add = Component(self.Xdata)
        self.add_component(comp_to_add,label=component_id)
        self._components[component_id] = comp_to_add
        #self._create_pixel_and_world_components(ndim=2)
        #self._shape = self.Xdata.shape
        
        
        #if self.hub and not is_present:
        #    msg = DataAddComponentMessage(self, component_id)
        #    self.hub.broadcast(msg)
        #    msg = ComponentsChangedMessage(self)
        #    self.hub.broadcast(msg)
            
        # To avoid circular references when saving objects with references to
        # the data, we make sure that all Data objects have a UUID that can
        # uniquely identify them.
        self.uuid = str(uuid.uuid4())
              
    #@property
    #def label(self):
    #    return self._label
    
    #@label.setter
    #def label(self, value):
    #    if getattr(self, '_label', None) != value:
    #        self._label = value
    #        self.broadcast(attribute='label')
    #    elif value is None:
    #        self._label = value
    
    
    #@property
    #def shape(self):
    #    return self.Xdata.shape
    
    #@property
    #def main_components(self):
    #    return self._main_components
        
    def get_kind(self, cid):
        return 'numerical' #FIX -- we should not assume X is the only cid?
        
    #def get_data(self, cid, view=None):
    #    if view is not None:
    #        subset = self.Xdata[view] #Maybe this makes a copy without file backing?
    #    else:
    #        subset = self.Xdata[:,:]
    #    return subset.X #FIX -- we should not assume X is the only cid? 
    
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
        
    def compute_histogram(self, cids, weights=None, range=None, bins=None, log=None, subset_state=None):
        """
        Compute an n-dimensional histogram with regularly spaced bins.
        
        Currently this only implements 1-D histograms.
        
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
        """
        if len(cids) > 2:
            raise NotImplementedError()
        
        ndim = len(cids)
        
        x = self.get_data(cids[0])
        if isinstance(x, categorical_ndarray):
            x = x.codes
        if ndim > 1:
            y = self.get_data(cids[1])
            if isinstance(y, categorical_ndarray):
                y = y.codes
        
        if weights is not None:
            w = self.get_data(weights)
            if isinstance(w, categorical_ndarray):
                w = w.codes
        else:
            w = None
        
        if subset_state is not None:
            mask = self.get_mask(subset_state)#.to_mask(self)
            x = x[mask]
            if ndim > 1:
                y = y[mask]
            if w is not None:
                w = w[mask]
        
        if ndim == 1:
            xmin, xmax = range[0]
            xmin, xmax = sorted((xmin, xmax))
            keep = (x >= xmin) & (x <= xmax)
        else:
            (xmin, xmax), (ymin, ymax) = range
            xmin, xmax = sorted((xmin, xmax))
            ymin, ymax = sorted((ymin, ymax))
            keep = (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)
        
        x = x[keep]
        if ndim > 1:
            y = y[keep]
        if w is not None:
            print(keep)
            w = w[keep]
        
        if len(x) == 0:
            return np.zeros(bins)
        
        if ndim > 1 and len(y) == 0:
            return np.zeros(bins)
        
        if log is not None and log[0]:
            if xmin < 0 or xmax < 0:
                return np.zeros(bins)
            xmin = np.log10(xmin)
            xmax = np.log10(xmax)
            x = np.log10(x)
        
        if ndim > 1 and log is not None and log[1]:
            if ymin < 0 or ymax < 0:
                return np.zeros(bins)
            ymin = np.log10(ymin)
            ymax = np.log10(ymax)
            y = np.log10(y)
        
        # By default fast-histogram drops values that are exactly xmax, so we
        # increase xmax very slightly to make sure that this doesn't happen, to
        # be consistent with np.histogram.
        if ndim >= 1:
            xmax += 10 * np.spacing(xmax)
        if ndim >= 2:
            ymax += 10 * np.spacing(ymax)
        
        if ndim == 1:
            range = (xmin, xmax)
            return histogram1d(x, range=range, bins=bins[0], weights=w)
        elif ndim > 1:
            range = [(xmin, xmax), (ymin, ymax)]
            return histogram2d(x, y, range=range, bins=bins, weights=w)