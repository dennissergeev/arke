# -*- coding: utf-8 -*-
"""
I/O functions
"""
import iris

from .coords import REDUNDANT_COORDS


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
