# -*- coding: utf-8 -*-
import copy
from datetime import datetime
import glob
import numpy as np
import os
import scipy.interpolate as scinter

import xray

from . import utils

def load_data(filenames, replace_unknowns=False, default_path='', verbose=False):
    if hasattr(filenames, '__iter__'):
        filenames = tuple([os.path.join(default_path, ifilename) for ifilename in filenames])
    else:
        filenames = glob.glob(os.path.join(default_path, filenames))

    datasets = []
    for fname in filenames:
        if verbose:
            print('Reading '+fname)
            print()
                
        dataset = xray.open_dataset(fname)
        if replace_unknowns:
            replace_unknown_names(dataset)
        if verbose:
            print(dataset)
            print()
            print()
        datasets.append(dataset)
            
    if hasattr(filenames, '__iter__'):
        return datasets
    else:
        return datasets[0]
    
    
def replace_unknown_names(dataset, default_name='unknown'):
    """ Replace missing names within a dataset by using STASH attribute """
    for ivar in dataset:
        newnames_dict = dict()
        if default_name in dataset[ivar].name:
            try:
                stash_id = dataset[ivar].attrs['um_stash_source']
                newname = stash_id_to_name(stash_id)
                newnames_dict[dataset[ivar].name] = newname
            except KeyError:
                print('Unable to rename '+dataset[ivar].name+': STASH attribute is missing; skipping')
                
        dataset.rename(newnames_dict, inplace=True)


def get_model_real_coords(vrbl, dims='tzyx'):
    pref_coords = dict(x=('x', 'grid_longitude', 'longitude'),
                       y=('y', 'grid_latitude', 'latitude'), 
                       z=('height', 'level_height', 'pressure', 'atmosphere_hybrid_height_coordinate'),
                       t=('time'))
    model_coords = []
    for iax in dims:
        for icoord in vrbl.coords:
            if icoord in pref_coords[iax]:
                model_coords.append(vrbl.coords[icoord])
                
    if len(vrbl.shape) != len(model_coords):
        print('WARNING! Number of coordinates does not match the input variable shape!')
    return model_coords
        

def regrid_model_to_obs(datacontainer, obs_coord, model_coords=None, obs_time_convert=True, rot_ll=True, shift=None, dims='tzyx'):
    if isinstance(obs_coord, tuple):
        if len(dims) != len(obs_coord):
            raise ValueError('Shape of the obs_coord does not equal to the dims keyword')
        else:
            obs_coord_dict = dict()
            for iax, i_coord in zip(dims, obs_coord):
                obs_coord_dict[iax] = i_coord
    elif isinstance(obs_coord, dict):
        obs_coord_dict = copy.deepcopy(obs_coord)
    else:
        raise ValueError('obs_coord can only be a tuple or a dict')    
    
    #print(np.nanmin(obs_coord_dict['x']),np.nanmax(obs_coord_dict['x']))
    #print(np.nanmin(obs_coord_dict['y']),np.nanmax(obs_coord_dict['y']))
    #print()
    #[ii+jj for ii, jj in zip(obs_coord[1:], zyxsh)]
    
    if isinstance(shift, dict):
        for iax in shift:
            obs_coord_dict[iax] += shift[iax]   
        
    ivar = dataset[varname]
    if model_coords is None:
        model_coords = get_model_real_coords(ivar, dims=dims)
    model_coord_points = [i.values for i in model_coords]
    
    if obs_time_convert:
        # Model times: converting numpy.datetime64 type array to datetimes,
        #              subtracting reference timestamp, converting to total number of seconds
        # Obs. times: subtracting ref. timestamp and converting to total seconds too
        #
        model_time = model_coords[0]
        ref_datetime = model_time.forecast_reference_time.values.astype(np.dtype('<M8[us]')).astype(datetime)
        model_coord_points[0] = np.array([(i.values.astype(np.dtype('<M8[us]')).astype(datetime) - ref_datetime).total_seconds() for i in model_time])
        obs_coord_dict['t'] = np.array([(i - ref_datetime).total_seconds() for i in obs_coord_dict['t']])
        #obs_tcoord = []
        #for idt in obs_coord[0]:
        #    obs_tcoord.append(np.datetime64(idt))
        #obs_coord[0] = np.array(obs_tcoord)
        
    if rot_ll:
        nplon = dataset.rotated_latitude_longitude.grid_north_pole_longitude
        nplat = dataset.rotated_latitude_longitude.grid_north_pole_latitude
        obs_coord_dict['x'], obs_coord_dict['y'] = iris.analysis.cartography.rotate_pole(obs_coord_dict['x'], obs_coord_dict['y'],
                                                                                         nplon, nplat)
    if isinstance(ivar, xray.dataset):
        ivar_data = ivar.data
    elif isinstance(ivar, np.ndarray):
        ivar_data = ivar

    fill_value = np.array(np.nan).astype(ivar.dtype)

    InterpFun = scinter.RegularGridInterpolator(model_coord_points, ivar_data, bounds_error=False, fill_value=fill_value)
    
    obs_coord_interp_arg = []
    for iax in dims:
        obs_coord_interp_arg.append(obs_coord_dict[iax])
    obs_coord_interp_arg = np.vstack(obs_coord_interp_arg).T
    #print(model_coord_points)
    #print()
    #print(obs_coord_interp_arg)
    return InterpFun(obs_coord_interp_arg)
