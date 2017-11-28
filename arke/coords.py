# -*- coding: utf-8 -*-
"""
Cube coords and related functions
"""
import iris
from iris.util import broadcast_to_shape
from iris.fileformats.pp import EARTH_RADIUS

PHYS_COORDS = dict(x=('x', 'grid_longitude', 'longitude'),
                   y=('y', 'grid_latitude', 'latitude'),
                   z=('height', 'level_height', 'pressure',
                      'atmosphere_hybrid_height_coordinate'),
                   t=('time'))
REDUNDANT_COORDS = ('forecast_period', )


def get_cube_datetimes(cube):
    """Get a list of `iris.cube.Cube`'s time points as `datetime.datetime`s """
    return cube.coord('time').units.num2date(cube.coord('time').points)


def nearest_plevel(cube, pval):
    """
    Find the value nearest to the given `pval`
    among points of pressure coordinate of a given cube
    """
    pcoord = cube.coord('pressure')
    i = pcoord.nearest_neighbour_index(pval)
    return pcoord.points[i]


def nearest_tval(cube, dt):
    """
    Find the value nearest to the given `dt` (`datetime.datetime` obj)
    among points of time coordinate of a given cube
    """
    timevar = cube.coord('time')
    itime = timevar.nearest_neighbour_index(timevar.units.date2num(dt))
    return timevar.units.num2date(timevar.points[itime])


def add_coord_system(src_cube, cs=None, inplace=False):
    """Add coordinate system to x and y coordinates of a cube"""
    if cs is None:
        cs = iris.coord_systems.GeogCS(EARTH_RADIUS)
    if inplace:
        cube = src_cube
    else:
        cube = src_cube.copy()
    for ax in ('x', 'y'):
        icoord = cube.coord(axis=ax)
        icoord.coord_system = cs
        # idim = cube.coord_dims(icoord)
        # cube.remove_coord(icoord)
        # icoord.coord_system = cs
        # cube.add_dim_coord(icoord, idim)
    if not inplace:
        return cube


def get_model_real_coords(cube, dims='tzyx', subtract360=True):
    """
    Retrieve 'model' coordinates of a cube along given dimensions

    Parameters
    ----------
    cube: iris.cube.Cube
        iris cube with coordinates
    dims: str, optional
        dimensions (axes) to do the search for
    subtract360: bool, optional (default True)
        subtract 360 from horizontal coordinates if they are all > 180
    Returns
    -------
    list
        'model' coordinates
    """
    model_coords = []
    for iax in dims:
        idim = cube.coords(axis=iax)
        if len(idim) > 1:
            for icoord in idim:
                if icoord.name() in PHYS_COORDS[iax]:
                    if (iax in 'xy' and all(icoord.points > 180.0)
                       and subtract360):
                        icoord.points = icoord.points - 360.0
                    model_coords.append(icoord)
        elif len(idim) == 1:
            if iax in 'xy' and all(idim[0].points > 180.0) and subtract360:
                idim[0].points = idim[0].points - 360.0
            model_coords.append(idim[0])

    if len(cube.shape) != len(model_coords):
        print('WARNING!'
              'Number of coordinates does not match the input variable shape!')
    return model_coords


def get_phys_coord(cube, axis, subtract360=True):
    """
    Retrieve a 'physical' coordinate along the given axis
    Parameters
    ----------
    cube: iris.cube.Cube
        iris cube with coordinates
    axis: str, optional
        dimension to do the search for
    subtract360: bool, optional (default True)
        subtract 360 from horizontal coordinates if they are all > 180
    Returns
    -------
    iris.coord.DimCoord
        'physical' coordinate
    """
    result = None
    coords = cube.coords(axis=axis)
    for icoord in coords:
        if icoord.name() in PHYS_COORDS[axis]:
            if axis in 'xy' and all(icoord.points > 180.0) and subtract360:
                icoord.points = icoord.points - 360.0
            result = icoord
            break

    return result


def pres_coord_to_cube(other_cube):
    """
    Convert pressure coordinate points to a cube of the same dimension as
    the given cube
    """
    pcoord = other_cube.coord('pressure')
    dim_map = other_cube.coord_dims('pressure')
    p_data = pcoord.points
    if len(dim_map) > 0:
        p_data = broadcast_to_shape(p_data,
                                    other_cube.shape,
                                    dim_map)
        dc = [(c.copy(), other_cube.coord_dims(c))
              for c in other_cube.dim_coords]
        ac = [(c.copy(), other_cube.coord_dims(c))
              for c in other_cube.aux_coords]
        pcube = iris.cube.Cube(data=p_data,
                               units=pcoord.units,
                               standard_name='air_pressure',
                               dim_coords_and_dims=dc,
                               aux_coords_and_dims=ac,
                               )
    else:
        pcube = iris.cube.Cube(data=p_data,
                               standard_name='air_pressure',
                               units=pcoord.units)
    return pcube
