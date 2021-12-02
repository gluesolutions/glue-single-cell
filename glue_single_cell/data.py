"""
A glue dataclass that wraps an AnnData object

The primary motivation is to support on-disk access
to a sparse datatype

There is a lot in Data that will probably just work normally
here. Maybe we should sub-class that instead?

"""
import pandas as pd
import numpy as np

from glue.core.data import BaseCartesianData
from glue.core.component_id import ComponentID, ComponentIDDict

from fast_histogram import histogram1d, histogram2d

import anndata

class DataAnnData(BaseCartesianData):
    def __init__(self, data, name):
        """
        data is an AnnData object with a single sparse
        matrix and (generally) with a filename backing
        """
        super(DataAnnData, self).__init__()
        print(data)
        self.Xdata = data
        #self.Xdata = data.X
        self.name = name
        self.id = ComponentIDDict(self)
        self._main_components = [ComponentID(label='X',parent=self)] #FIX -- we should not assume X is the only cid?
        
        #Coords are like this:
        adata.obs_names.astype('int')
        
    @property
    def label(self):
        return self.name
    
    @property
    def shape(self):
        return self.Xdata.shape
    
    @property
    def main_components(self):
        return self._main_components
        
    def get_kind(self, cid):
        return 'numerical' #FIX -- we should not assume X is the only cid?
        
    def get_data(self, cid, view=None):
        if view is not None:
            subset = self.Xdata[view] #Maybe this makes a copy without file backing?
        else:
            subset = self.Xdata[:,:]
        return subset.X #FIX -- we should not assume X is the only cid? 
    
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
        
        if ndim == 1:
                xmin, xmax = range[0]
                xmin, xmax = sorted((xmin, xmax))
                keep = (x >= xmin) & (x <= xmax)
        else:
            (xmin, xmax), (ymin, ymax) = range
            xmin, xmax = sorted((xmin, xmax))
            ymin, ymax = sorted((ymin, ymax))
            keep = (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)
        
        if x.dtype.kind == 'M':
            x = datetime64_to_mpl(x)
            xmin = datetime64_to_mpl(xmin)
            xmax = datetime64_to_mpl(xmax)
        else:
            keep &= ~np.isnan(x)
        
        if ndim > 1:
            if y.dtype.kind == 'M':
                y = datetime64_to_mpl(y)
                ymin = datetime64_to_mpl(ymin)
                ymax = datetime64_to_mpl(ymax)
            else:
                keep &= ~np.isnan(y)
        print("Got to this point in the middle of compute_histogram")
        x = x[keep]
        if ndim > 1:
            y = y[keep]
        if w is not None:
            w = w[keep]
            
        
        if ndim == 1:
            range = (xmin, xmax)
            return histogram1d(x, range=range, bins=bins[0], weights=w)
        elif ndim > 1:
            range = [(xmin, xmax), (ymin, ymax)]
            print("Trying to compute histogram2d")
            print(f"x = {x}")
            print(f"y = {y}")
            print(f"range={range}")
            return histogram2d(x, y, range=range, bins=bins, weights=w)


