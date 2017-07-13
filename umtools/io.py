# -*- coding: utf-8 -*-
"""
I/O functions
"""
import iris

from .coords import REDUNDANT_COORDS
from .numerics import (nearest_xy_grid_2d_index,
                       mask_cube_outside_circle_xy)


def slice_cubelist(cubelist, coord_name, coord_values):
    """
    Filter cubelist by coordinate name and values
    iterating over cubes and using the nearest coordinate value.

    If a coordinate is not found, include in the returned cubelist anyway.
    """
    if isinstance(coord_values, (int, float)):
        coord_values = [coord_values]
    extracted_cubelist = iris.cube.CubeList()
    for cube in cubelist:
        try:
            coord = cube.coord(coord_name)
            idx = []
            for coord_value in coord_values:
                idx.append(coord.nearest_neighbour_index(coord_value))
            constr = {coord.name(): coord.points[idx]}
            extracted_cubelist.append(cube.extract(iris.Constraint(**constr)))
        except iris.exceptions.CoordinateNotFoundError:
            extracted_cubelist.append(cube)
    return extracted_cubelist


def extract_levels(cubelist, level_dict):
    vert_subset = dict()
    for zname in set([d['name'] for d in level_dict.values() if d['name']]):
        vert_subset[zname] = []
    for v in level_dict.values():
        if v['name'] in vert_subset.keys():
            vert_subset[v['name']].append(v['value'])

    extracted = cubelist
    for coord_name, coord_values in vert_subset.items():
        extracted = slice_cubelist(extracted, coord_name, coord_values)
    return extracted


def subset_cubelist(cubelist, h_subset):
    cl = iris.cube.CubeList()
    for cube in cubelist:
        if 'ilon' in h_subset and 'ilat' in h_subset:
            iy, ix = nearest_xy_grid_2d_index(cube,
                                              h_subset['ilat'],
                                              h_subset['ilon'])
        elif 'ix' in h_subset and 'iy' in h_subset:
            iy, ix = h_subset['iy'], h_subset['ix']
        else:
            iy, ix = [i//2 for i in (cube.shape[-2], cube.shape[-1])]
        # print('subset_cubelist:', cube.name())
        # print(cube.shape)
        # print(iy, ix)

        if h_subset['method'] == 'wh':
            xslice = slice(ix-h_subset['w']//2,
                           ix+h_subset['w']//2+1)
            yslice = slice(iy-h_subset['y']//2,
                           iy+h_subset['y']//2+1)
            cl.append(cube[..., yslice, xslice])
        elif h_subset['method'] == 'xy':
            xslice = slice(h_subset['corners'][0], h_subset['corners'][1])
            yslice = slice(h_subset['corners'][2], h_subset['corners'][3])
            cl.append(cube[..., yslice, xslice])
        elif h_subset['method'] == 'radius':
            # print(cube, h_subset['r'], ix, iy)
            dx = cube.attributes['um_res'].to_flt('km')
            mc = mask_cube_outside_circle_xy(cube, h_subset['r'], ix, iy,
                                             dx=dx,
                                             return_copy=True,
                                             return_mask=False)
            # xslice = yslice = slice(None)
            cl.append(mc)
        else:
            raise NotImplementedError()

    return cl


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


def clean_call(cube, field, filename):
    try:
        for factory in cube.aux_factories:
            cube.remove_aux_factory(factory)
    except:
        pass
    try:
        cube.remove_coord(cube.coord('altitude'))
    except:
        pass
    try:
        cube.remove_coord(cube.coord('surface_altitude'))
    except:
        pass
