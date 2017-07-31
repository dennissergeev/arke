# -*- coding: utf-8 -*-
"""
Regridding and interpolation
"""
import copy
import iris
import numpy as np
import scipy.interpolate as scinter

from . import io, coords


def interp_cubes_to_points(cubelist, cube_name_and_axes,
                           verbose=0, extrapolation_mode='linear'):
    """ """
    scheme = iris.analysis.Linear(extrapolation_mode=extrapolation_mode)
    src_cube = io.get_cube(cubelist, cube_name_and_axes['name'])
    pnts = []
    if 'dim_ave' in cube_name_and_axes:
        for iax in cube_name_and_axes['dim_ave']:
            iax_name = src_cube.coord(axis=iax).name()
            iax_pnts = 0.5 * (src_cube.coord(axis=iax).points[1:] +
                              src_cube.coord(axis=iax).points[:-1])
            if all(iax_pnts > 180.0):
                iax_pnts = iax_pnts - 360.0
            pnts.append((iax_name, iax_pnts))

    else:
        for iax in 'xy':
            iax_name = src_cube.coord(axis=iax).name()
            iax_pnts = src_cube.coord(axis=iax).points
            if all(iax_pnts > 180.0):
                iax_pnts = iax_pnts - 360.0
            pnts.append((iax_name, iax_pnts))

    for k, icube in enumerate(cubelist):
        new_cube = icube.interpolate(pnts, scheme=scheme)
        if verbose > 1:
            print('Interpolation of {} is completed.'.format(new_cube.name()))
        cubelist[k] = new_cube


def regrid_model_to_obs(datacontainer, obs_coord, model_coords=None,
                        obs_time_convert=True, rot_ll=True, shift=None,
                        dims='tzyx', return_cube=False):
    if isinstance(obs_coord, tuple):
        if len(dims) != len(obs_coord):
            _msg = 'Shape of the obs_coord does not equal to the dims keyword'
            raise ValueError(_msg)
        else:
            obs_coord_dict = dict()
            for iax, i_coord in zip(dims, obs_coord):
                obs_coord_dict[iax] = i_coord
    elif isinstance(obs_coord, dict):
        obs_coord_dict = copy.deepcopy(obs_coord)
    else:
        raise ValueError('obs_coord can only be a tuple or a dict')

    # print(np.nanmin(obs_coord_dict['x']),np.nanmax(obs_coord_dict['x']))
    # print(np.nanmin(obs_coord_dict['y']),np.nanmax(obs_coord_dict['y']))
    # print()
    # [ii+jj for ii, jj in zip(obs_coord[1:], zyxsh)]

    if isinstance(shift, dict):
        for iax in shift:
            obs_coord_dict[iax] += shift[iax]

    if isinstance(datacontainer, (list, tuple)):
        ivar = io.get_cube(datacontainer[0], datacontainer[1])
    elif isinstance(datacontainer, dict):
        ivar = io.get_cube(datacontainer['dataset'], datacontainer['varname'])
    elif isinstance(datacontainer, iris.cube.Cube):
        ivar = datacontainer
    elif isinstance(datacontainer, np.ndarray):
        if model_coords is None or not rot_ll:
            _msg = 'Model coords must be passed explicitly'\
                   ' if the input data is numpy.ndarray'
            raise ValueError(_msg)
        else:
            ivar = datacontainer
    else:
        raise ValueError('Unrecognized input data type')
    if model_coords is None:
        model_coords = coords.get_model_real_coords(ivar, dims=dims)
    model_coord_points = [i.points for i in model_coords]

    if obs_time_convert and 't' in dims:
        obs_coord_dict['t'] = (model_coords[0].
                               units.
                               date2num(obs_coord_dict['t']))

    if rot_ll:
        try:
            nplon = ivar.coord_system().grid_north_pole_longitude
            nplat = ivar.coord_system().grid_north_pole_latitude
        except:
            nplon = model_coords[-1].coord_system.grid_north_pole_longitude
            nplat = model_coords[-1].coord_system.grid_north_pole_latitude
        obs_coord_dict['x'],
        obs_coord_dict['y'] = (iris.analysis.cartography.
                               rotate_pole(obs_coord_dict['x'],
                                           obs_coord_dict['y'],
                                           nplon, nplat))
    # print(np.nanmin(obs_coord_dict['x']),np.nanmax(obs_coord_dict['x']))
    # print(np.nanmin(obs_coord_dict['y']),np.nanmax(obs_coord_dict['y']))
    if isinstance(ivar, iris.cube.Cube):
        ivar_data = ivar.data
    elif isinstance(ivar, np.ndarray):
        ivar_data = ivar

    fill_value = np.array(np.nan).astype(ivar.dtype)

    InterpFun = scinter.RegularGridInterpolator(model_coord_points, ivar_data,
                                                bounds_error=False,
                                                fill_value=fill_value)

    obs_coord_interp_arg = []
    for iax in dims:
        obs_coord_interp_arg.append(obs_coord_dict[iax])
    obs_coord_interp_arg = np.vstack(obs_coord_interp_arg).T
    # print('-------------------------------------')
    # print(ivar.coord('time').units.num2date(model_coord_points[0]))
    # print()
    # print(ivar.coord('time').units.num2date(obs_coord_interp_arg[:, 0]))
    # print('-------------------------------------')
    interp_array = InterpFun(obs_coord_interp_arg)

    if return_cube:
        interp_cube = iris.cube.Cube(interp_array)
        interp_cube.rename(u'interpolated_'+ivar.name())
        interp_cube.units = ivar.units
        if 't' in dims:
            numeric_time = obs_coord_dict['t']
            t = iris.coords.DimCoord(numeric_time, 'time',
                                     units=ivar.coord('time').units)
            interp_cube.add_dim_coord(t, 0)
        if 'x' in dims:
            x = iris.coords.AuxCoord(obs_coord['x'], 'longitude')
            interp_cube.add_aux_coord(x, 0)
        if 'y' in dims:
            y = iris.coords.AuxCoord(obs_coord['y'], 'latitude')
            interp_cube.add_aux_coord(y, 0)
        if 'z' in dims:
            z = iris.coords.AuxCoord(obs_coord['z'], 'height')
            interp_cube.add_aux_coord(z, 0)
        return interp_cube
    else:
        return interp_array
