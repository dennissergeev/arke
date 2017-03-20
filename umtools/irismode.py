# -*- coding: utf-8 -*-
import configparser
import copy
import glob
import numpy as np
import os
import scipy.interpolate as scinter

import cf_units
import iris
from iris.fileformats.pp import EARTH_RADIUS

from . import utils

PHYS_COORDS = dict(x=('x', 'grid_longitude', 'longitude'),
                   y=('y', 'grid_latitude', 'latitude'),
                   z=('height', 'level_height', 'pressure',
                      'atmosphere_hybrid_height_coordinate'),
                   t=('time'))
REDUNDANT_COORDS = ('forecast_period', )


def add_coord_system(src_cube, cs=None):
    if cs is None:
        cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
    cube = src_cube.copy()
    for ax in ('x', 'y'):
        icoord = cube.coord(axis=ax)
        idim = cube.coord_dims(icoord)
        cube.remove_coord(icoord)
        icoord.coord_system = cs
        cube.add_dim_coord(icoord, idim)
    return cube


def interp_cubes_to_points(cubelist, cube_name_and_axes,
                           verbose=0, extrapolation_mode='linear'):
    """ """
    scheme = iris.analysis.Linear(extrapolation_mode=extrapolation_mode)
    src_cube = get_cube(cubelist, cube_name_and_axes['name'])
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
        if verbose > 1:
            utils.tic()
        new_cube = icube.interpolate(pnts, scheme=scheme)
        if verbose > 1:
            print('Interpolation of {} is completed.'.format(new_cube.name()))
            utils.toc()
            print()
        cubelist[k] = new_cube


def load_data(filenames,
              replace_unknowns=False, default_path='', verbose=False):
    if isinstance(filenames, list) or isinstance(filenames, tuple):
        filenames = tuple([os.path.join(default_path, ifilename)
                           for ifilename in filenames])
    elif isinstance(filenames, str):
        filenames = sorted(glob.glob(os.path.join(default_path, filenames)))
    else:
        raise TypeError('filenames should be a string or a list of strings')

    if verbose:
        print('File names list: {}'.format(filenames))
    datasets = []
    for fname in filenames:
        if verbose:
            print('Reading '+fname)
            print()

        cubelist = iris.load(fname)
        if replace_unknowns:
            replace_unknown_names(cubelist)
        if verbose:
            print(cubelist)
            print()
            print()
        datasets.append(cubelist)

    if isinstance(filenames, list) or isinstance(filenames, tuple):
        return datasets
    else:
        return datasets[0]


def replace_unknown_names(dataset, default_name='unknown'):
    """ Replace missing names within an `iris.cube.CubeList`
        by using STASH attribute """
    for ivar in dataset:
        if default_name in ivar.name().lower():
            try:
                stash_id = ivar.attributes['STASH'].__str__()
                ivar.rename(utils.stash_id_to_name(stash_id))
            except AttributeError:
                print('Unable to rename, STASH attribute is missing')


def get_cube(cubelist, cube_name, lazy=True):
    """ Return a `cube` from a `iris.cube.CubeList` by name` """
    i = None
    for i in cubelist:
        if lazy:
            match = cube_name.lower() in i.name().lower()
        else:
            match = cube_name == i.name()

        if match:
            return i
    if i is None:
        _msg = 'Cube with name {0} not found in {1}'
        raise ValueError(_msg.format(cube_name, cubelist))


def get_model_real_coords(vrbl, dims='tzyx'):
    """ Retrieve 'physical' coordinates of """
    model_coords = []
    for iax in dims:
        idim = vrbl.coords(axis=iax)
        if len(idim) > 1:
            for icoord in idim:
                if icoord.name() in PHYS_COORDS[iax]:
                    if iax in 'xy' and all(icoord.points > 180.0):
                        icoord.points = icoord.points - 360.0
                    model_coords.append(icoord)
        elif len(idim) == 1:
            if iax in 'xy' and all(idim[0].points > 180.0):
                idim[0].points = idim[0].points - 360.0
            model_coords.append(idim[0])

    if len(vrbl.shape) != len(model_coords):
        print('WARNING!'
              'Number of coordinates does not match the input variable shape!')
    return model_coords


def get_phys_coord(vrbl, axis, subtract360=True):
    """ Retrieve 'physical' coordinates of """
    result = None
    coords = vrbl.coords(axis=axis)
    for icoord in coords:
        if icoord.name() in PHYS_COORDS[axis]:
            if axis in 'xy' and all(icoord.points > 180.0) and subtract360:
                icoord.points = icoord.points - 360.0
            result = icoord
            break

    return result


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
        ivar = get_cube(datacontainer[0], datacontainer[1])
    elif isinstance(datacontainer, dict):
        ivar = get_cube(datacontainer['dataset'], datacontainer['varname'])
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
        model_coords = get_model_real_coords(ivar, dims=dims)
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


def unrotate_uv(u, v, target_cs=None, remove_aux_xy=True):
    if target_cs is None:
        target_cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
    uv = iris.analysis.cartography.rotate_winds(u, v, target_cs)
    if remove_aux_xy:
        [cube.remove_coord(i) for i in ('projection_x_coordinate',
                                        'projection_y_coordinate')
         for cube in uv]
    return uv


def unrotate_wind(cubelist,
                  uwind_name='x_wind', vwind_name='y_wind',
                  newcs=iris.coord_systems.GeogCS(EARTH_RADIUS),
                  replace=False, verbose=0):

        u = get_cube(cubelist, uwind_name, lazy=False)
        v = get_cube(cubelist, vwind_name, lazy=False)

        if u is not None or v is not None:
            oldcs = u.coord_system()
            if verbose > 1:
                print('Rotating winds from {} to {}'.format(oldcs, newcs))
                print()
            u_rot, v_rot = iris.analysis.cartography.rotate_winds(u, v, newcs)
            if replace:
                cubelist[cubelist.index(u)] = u_rot
                cubelist[cubelist.index(v)] = v_rot
            else:
                cubelist.append(u_rot)
                cubelist.append(v_rot)
        else:

            if get_cube(cubelist, 'transformed_x_wind',
                        lazy=False) is not None\
                and get_cube(cubelist, 'transformed_y_wind',
                             lazy=False) is not None:
                print('transformed winds are in the file already')
            else:
                print('u-wind or v-wind cubes not found. No winds rotating.')


def convert_unit_str(str1, str2):
    return cf_units.Unit(str1).convert(1, cf_units.Unit(str2))


def extract_as_single_cube(cubelist, constraints):
    try:
        cube = cubelist.extract(constraints, strict=True)
    except iris.exceptions.ConstraintMismatchError:
        cube = None
        cubes = cubelist.extract(constraints)
        for icoord in REDUNDANT_COORDS:
            conc = cubes.concatenate()
            if len(conc) == 1:
                cube = conc[0]
                break
            else:
                for icube in cubes:
                    icube.remove_coord(icoord)
                conc = cubes.concatenate()
                if len(conc) == 1:
                    cube = conc[0]
                    break
        if cube is None:
            raise ValueError('Unable to concatenate')
    return cube


def get_cube_datetimes(cube):
    return cube.coord('time').units.num2date(cube.coord('time').points)


def rename_cubes_from_stashmaster(cube, stashmaster_file,
                                  parser=configparser.ConfigParser(),
                                  add_help=True):
    if cube.name().lower() == 'unknown':
        stash = cube.attributes['STASH']
        parser.read(stashmaster_file)
        try:
            section = parser['stashmaster:code({})'.format(stash.item)]
            new_name = '_'.join(section['description'].split())
            cube.rename(new_name)
            if add_help:
                try:
                    cube.attributes['help'] = section['help'].\
                                              replace('=', '').\
                                              replace('\n', '')
                except KeyError:
                    pass
        except KeyError:
            pass
