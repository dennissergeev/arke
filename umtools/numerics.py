# -*- coding: utf-8 -*-
import functools
import iris
from iris.fileformats.pp import EARTH_RADIUS
import numpy as np
import operator

from . import io


def nearest_xy_grid_2d_index(cube, point_lon, point_lat):
    lons, lats = unrotate_xy_grids(cube)
    dist = np.sqrt((lons - point_lon)**2 + (lats - point_lat)**2)
    iy, ix = np.unravel_index(np.argmin(dist), lons.shape)

    return (iy, ix)


def mask_cube_outside_circle_xy(cube, radius, point_lon, point_lat,
                                dx=1, dy=None,
                                return_copy=True, return_mask=True):
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
        remainder = functools.reduce(operator.mul, shp)
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

    Args:

        * cube - The cube with rotated coordinate system
                 for which to generate 2D X and Y unrotated coordinates.

    Example::

        lon, lat = unrotate_xy_grids(cube)

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

        u = io.get_cube(cubelist, uwind_name, lazy=False)
        v = io.get_cube(cubelist, vwind_name, lazy=False)

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

            if io.get_cube(cubelist, 'transformed_x_wind',
                           lazy=False) is not None\
                and io.get_cube(cubelist, 'transformed_y_wind',
                                lazy=False) is not None:
                print('transformed winds are in the file already')
            else:
                print('u-wind or v-wind cubes not found. No winds rotating.')
