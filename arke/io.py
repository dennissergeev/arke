# -*- coding: utf-8 -*-
"""
I/O functions
"""
import iris

from .coords import REDUNDANT_COORDS


def extract_as_single_cube(cubelist, constraints):
    """
    Extract cube from a cubelist as a single cube

    If concatenation does not produce a single cube, try removing
    'redundant' coordinates first

    Parameters
    ----------
    cubelist: iris.cube.CubeList
        List of cubes
    constraints: iris.Constraint
        Constraint(s) used to extract cubes
    Returns
    -------
    iris.cube.Cube
        a single cube
    """
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
            _msg = 'Unable to concatenate {}'.format(cubes)
            raise iris.exceptions.ConcatenateError([_msg])
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
    """Clean some of the cube coordinates when loading data"""
    try:
        for factory in cube.aux_factories:
            cube.remove_aux_factory(factory)
    except iris.exceptions.CoordinateNotFoundError:
        pass
    try:
        cube.remove_coord(cube.coord('altitude'))
    except iris.exceptions.CoordinateNotFoundError:
        pass
    try:
        cube.remove_coord(cube.coord('surface_altitude'))
    except iris.exceptions.CoordinateNotFoundError:
        pass
