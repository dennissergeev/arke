# -*- coding: utf-8 -*-
"""
Regridding and interpolation
"""
import copy
from functools import reduce
import iris
from iris.fileformats.pp import EARTH_RADIUS
from iris.analysis.cartography import rotate_winds
import numpy as np
import operator
import scipy.interpolate as scinter

from . import io, coords
from .exceptions import ArgumentError, CoordinateRankError


def nearest_xy_grid_2d_index(cube, point_lon, point_lat):
    """
    Find index-pair nearest to the given longitude an latitude

    Parameters
    ----------
    cube: iris.cube.Cube
        iris cube
    point_lon: float
        longitude of the point
    point_lat: float
        latitude of the point
    Returns
    -------
    (int, int)
        a pair of indices
    """
    lons, lats = unrotate_xy_grids(cube)
    dist = np.sqrt((lons - point_lon)**2 + (lats - point_lat)**2)
    iy, ix = np.unravel_index(np.argmin(dist), lons.shape)

    return (iy, ix)


def mask_cube_outside_circle_xy(cube, radius, point_lon, point_lat,
                                dx=1, dy=None,
                                return_copy=True, return_mask=True):
    """
    Mask a cube outside a circle

    Parameters
    ----------
    cube: iris.cube.Cube
        iris cube
    radius: float
        radius around `point_lon` and `point_lat`
    point_lon: float
        longitude of the point
    point_lat: float
        latitude of the point
    dx: float, optional (default 1)
        grid spacing of the `cube`'s x-axis
        should be in the same units as radius
    dy: float, optional (default =dx)
        grid spacing in along y-axis
    return_copy: bool (default True)
        if True, return a masked copy of the given cube,
        otherwise mask the given cube
    return_mask: bool (default True)
        if True, return the mask too
    Returns
    -------
    iris.cube.Cube
        (if return_copy=True) masked cube
    """
    dy = dx if dy is None else dy
    iy, ix = nearest_xy_grid_2d_index(cube, point_lon, point_lat)
    xcoord, _ = iris.analysis.cartography.get_xy_grids(cube)
    ny, nx = xcoord.shape
    x, y = np.ogrid[-iy:ny-iy, -ix:nx-ix]
    xdist, ydist = x*dx, y*dy
    mask = xdist*xdist + ydist*ydist > radius*radius
    shp = cube.shape
    if mask.shape != shp:
        shp = list(shp)
        [shp.remove(n) for n in (nx, ny)]
        remainder = reduce(operator.mul, shp)
        mask = mask.repeat(remainder).reshape((ny, nx, *shp))
        mask = np.rollaxis(np.rollaxis(mask, 0, mask.ndim), 0, mask.ndim)
    masked_data = np.ma.masked_where(mask, cube.data)

    if return_copy:
        cp_cube = cube.copy()
        cp_cube.data = masked_data
        if return_mask:
            return cp_cube, mask
        else:
            return cp_cube
    else:
        cube.data = masked_data
        if return_mask:
            return mask


def unrotate_xy_grids(cube):
    """
    Convert rotated-pole lons and lats to unrotated ones
    using X and Y coordinate for a given cube.
    Unites get_xy_grids() and unrotate_pole() functions.

    Parameters
    ----------
    cube: iris.cube.Cube
        The cube with rotated coordinate system
        for which to generate 2D X and Y unrotated coordinates.
    Returns
    -------
    tuple of 2 numpy arrays (x, y)

    Example
    -------
    x, y = unrotate_xy_grids(cube)
    """
    x, y = iris.analysis.cartography.get_xy_grids(cube)

    cs = cube.coord_system('CoordSystem')
    if isinstance(cs, iris.coord_systems.RotatedGeogCS):
        nplon = cube.coord_system().grid_north_pole_longitude
        nplat = cube.coord_system().grid_north_pole_latitude
        x, y = iris.analysis.cartography.unrotate_pole(x, y, nplon, nplat)
    else:
        # no rotation needed
        pass

    return (x, y)


def unrotate_lonlat_grids(cube):
    """
    Convert rotated-pole lons and lats to unrotated ones
    using lon and lat coordinate for a given cube.
    Based on _get_lat_lon_coords() and unrotate_pole() functions.

    Parameters
    ----------
    cube: iris.cube.Cube
        The cube with rotated coordinate system
        for which to generate 2D X and Y unrotated coordinates.
    Returns
    -------
    tuple of 2 numpy arrays (lons, lats)

    Example
    -------
    lon, lat = unrotate_latlon_grids(cube)
    """
    y_coord, x_coord = iris.analysis.cartography._get_lat_lon_coords(cube)
    x = x_coord.points
    y = y_coord.points

    if x.ndim == y.ndim == 1:
        # Convert to 2D.
        x, y = np.meshgrid(x, y)
    elif x.ndim == y.ndim == 2:
        # They are already in the correct shape.
        pass
    else:
        raise CoordinateRankError("Expected 1D or 2D XY coords")

    cs = cube.coord_system('CoordSystem')
    if isinstance(cs, iris.coord_systems.RotatedGeogCS):
        nplon = cube.coord_system().grid_north_pole_longitude
        nplat = cube.coord_system().grid_north_pole_latitude
        x, y = iris.analysis.cartography.unrotate_pole(x, y, nplon, nplat)
    else:
        # no rotation needed
        pass

    return (x, y)


def unrotate_uv(u, v, target_cs=None, remove_aux_xy=True):
    """
    Rotate u- and v-winds to a given CS (Geog by default) and remove
    auxiliary coordinates created automatically by `rotate_winds()` function
    """
    if target_cs is None:
        target_cs = iris.coord_systems.GeogCS(EARTH_RADIUS)
    uv = rotate_winds(u, v, target_cs)
    if remove_aux_xy:
        [cube.remove_coord(i) for i in ('projection_x_coordinate',
                                        'projection_y_coordinate')
         for cube in uv]
    return uv


# def unrotate_wind(cubelist,
#                   uwind_name='x_wind', vwind_name='y_wind',
#                   newcs=iris.coord_systems.GeogCS(EARTH_RADIUS),
#                   replace=False, verbose=0):
#
#         u = io.get_cube(cubelist, uwind_name, lazy=False)
#         v = io.get_cube(cubelist, vwind_name, lazy=False)
#
#         if u is not None or v is not None:
#             oldcs = u.coord_system()
#             if verbose > 1:
#                 print('Rotating winds from {} to {}'.format(oldcs, newcs))
#                 print()
#             u_rot, v_rot = rotate_winds(u, v, newcs)
#             if replace:
#                 cubelist[cubelist.index(u)] = u_rot
#                 cubelist[cubelist.index(v)] = v_rot
#             else:
#                 cubelist.append(u_rot)
#                 cubelist.append(v_rot)
#         else:
#
#             if io.get_cube(cubelist, 'transformed_x_wind',
#                            lazy=False) is not None\
#                 and io.get_cube(cubelist, 'transformed_y_wind',
#                                 lazy=False) is not None:
#                 print('transformed winds are in the file already')
#             else:
#                 print('u-wind or v-wind cubes not found. No winds rotating.')


# def interp_cubes_to_points(cubelist, cube_name_and_axes,
#                            verbose=0, extrapolation_mode='linear'):
#     """Interpolate cube to  """
#     scheme = iris.analysis.Linear(extrapolation_mode=extrapolation_mode)
#     src_cube = io.get_cube(cubelist, cube_name_and_axes['name'])
#     pnts = []
#     if 'dim_ave' in cube_name_and_axes:
#         for iax in cube_name_and_axes['dim_ave']:
#             iax_name = src_cube.coord(axis=iax).name()
#             iax_pnts = 0.5 * (src_cube.coord(axis=iax).points[1:] +
#                               src_cube.coord(axis=iax).points[:-1])
#             if all(iax_pnts > 180.0):
#                 iax_pnts = iax_pnts - 360.0
#             pnts.append((iax_name, iax_pnts))
#
#     else:
#         for iax in 'xy':
#             iax_name = src_cube.coord(axis=iax).name()
#             iax_pnts = src_cube.coord(axis=iax).points
#             if all(iax_pnts > 180.0):
#                 iax_pnts = iax_pnts - 360.0
#             pnts.append((iax_name, iax_pnts))
#
#     for k, icube in enumerate(cubelist):
#         new_cube = icube.interpolate(pnts, scheme=scheme)
#         if verbose > 1:
#             print('Interpolation of {} completed.'.format(new_cube.name()))
#         cubelist[k] = new_cube


def regrid_model_to_obs(datacontainer, obs_coord, model_coords=None,
                        obs_time_convert=True, rot_ll=True, shift=None,
                        dims='tzyx', return_cube=False):
    # TODO: deprecate or rewrite the function
    if isinstance(obs_coord, tuple):
        if len(dims) != len(obs_coord):
            _msg = 'Shape of the obs_coord does not equal to the dims keyword'
            raise CoordinateRankError(_msg)
        else:
            obs_coord_dict = dict()
            for iax, i_coord in zip(dims, obs_coord):
                obs_coord_dict[iax] = i_coord
    elif isinstance(obs_coord, dict):
        obs_coord_dict = copy.deepcopy(obs_coord)
    else:
        raise ArgumentError('obs_coord can only be a tuple or a dict')

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
            raise ArgumentError(_msg)
        else:
            ivar = datacontainer
    else:
        raise ArgumentError('Unrecognized input data type')
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
        except AttributeError:
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
