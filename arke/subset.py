# -*- coding: utf-8 -*-
"""
Functions to subset cubes and cubelists
"""
import iris
from iris.analysis import trajectory
from iris.analysis.cartography import rotate_pole
import numpy as np

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


def extract_vert_section(cube, pnts):
    """
    Extract vertical slice of a cube

    Assume cube's dims: (..., y, x)
    """
    sample_count = pnts['sample_count']
    zcoord = cube.coord(axis='z', dim_coords=True)
    zspan = int(zcoord.points[0]), int(zcoord.points[-1])
    if 'x' in pnts and 'y' in pnts:
        # Grid indices (no interpolation needed)
        if len(pnts['x']) == 2 and len(pnts['y']) == 2:
            if ((pnts['x'] == pnts['x'][0]).all() and
               not (pnts['y'] == pnts['y'][0]).all()):
                # X const, Y varies
                xslice = pnts['x'][0]
                yslice = slice(*pnts['y'])
                sect_info = dict(v=zcoord.name(),
                                 h='y_axis_ind',
                                 label='z{}-{}y{}-{}x{}'.format(*zspan,
                                                                *pnts['y'],
                                                                pnts['x'][0]))
            elif ((pnts['y'] == pnts['y'][0]).all() and
                  not (pnts['x'] == pnts['x'][0]).all()):
                # Y const, X varies
                yslice = pnts['y'][0]
                xslice = slice(*pnts['x'])
                sect_info = dict(v=zcoord.name(),
                                 h='x_axis_ind',
                                 label='z{}-{}y{}x{}-{}'.format(*zspan,
                                                                pnts['y'][0],
                                                                *pnts['x']))
            sect = cube[..., yslice, xslice]
        else:
            _msg = (f'Only 2 point pairs is allowed for grid indices option'
                    ',\n{pnts} are given')
            raise NotImplementedError(_msg)
    elif 'lon' in pnts and 'lat' in pnts:
        # Lon / lat coordinates
        cs = cube.coord_system()
        xcoord = cube.coord(axis='x', dim_coords=True)
        ycoord = cube.coord(axis='y', dim_coords=True)
        xp, yp = xcoord.points, ycoord.points
        if (xp > 180).any():
            xp = xp - 360

        if isinstance(cs, (iris.coord_systems.RotatedGeogCS, )):
            # Rotated coordinate system
            # Use iris interpolation along a trajectory
            nplat = cs.grid_north_pole_latitude
            nplon = cs.grid_north_pole_longitude
            waypoints = []
            for ilon, ilat in zip(pnts['lon'], pnts['lat']):
                waypoints.append({'longitude': ilon, 'latitude': ilat})
            if sample_count == 'auto':
                _x, _y = [], []
                for ilon, ilat in zip(pnts['lon'], pnts['lat']):
                    _ilon, _ilat = rotate_pole(np.asarray(ilon),
                                               np.asarray(ilat),
                                               nplon, nplat)
                    _x.append(np.argmin(abs(xp - _ilon[0])))
                    _y.append(np.argmin(abs(yp - _ilat[0])))
                sample_count = int(((np.diff(_x)**2 +
                                     np.diff(_y)**2)**0.5).sum())
            else:
                assert isinstance(sample_count, int)

            traj = trajectory.Trajectory(waypoints, sample_count=sample_count)
            lon = np.array([d['longitude'] for d in traj.sampled_points])
            lat = np.array([d['latitude'] for d in traj.sampled_points])
            rlon, rlat = iris.analysis.cartography.rotate_pole(lon, lat,
                                                               nplon, nplat)
            sampled_points = [(xcoord.name(), rlon),
                              (ycoord.name(), rlat)]

            sect = trajectory.interpolate(cube, sampled_points)
            _add_real_lonlat(sect, lon, lat)

            latspan = pnts['lat'][0], pnts['lat'][-1]
            lonspan = pnts['lon'][0], pnts['lon'][-1]
            sect_info = dict(v=zcoord.name(),
                             h='index',
                             label='z{}-{}lat{}-{}lon{}-{}'.format(*zspan,
                                                                   *latspan,
                                                                   *lonspan))

        elif isinstance(cs, (None, iris.coord_systems.GeogCS)):
            # Rectangular grid
            if len(pnts['lon']) == 2 and len(pnts['lat']) == 2:
                if ((pnts['lon'] == pnts['lon'][0]).all() and
                   not (pnts['lat'] == pnts['lat'][0]).all()):
                    # X const, Y varies
                    xslice = np.argmin(abs(xp - pnts['lon'][0]))
                    yslice = slice(*[np.argmin(abs(yp - i))
                                     for i in pnts['lat']])

                elif ((pnts['lat'] == pnts['lat'][0]).all() and
                      not (pnts['lon'] == pnts['lon'][0]).all()):
                    # Y const, X varies
                    yslice = np.argmin(abs(yp - pnts['lat'][0]))
                    xslice = slice(*[np.argmin(abs(xp - i))
                                     for i in pnts['lon']])
                sect = cube[..., yslice, xslice]
                sect_info = {}  # TODO
            else:
                # Diagonal
                raise NotImplementedError()
        else:
            raise NotImplementedError(f"Can't deal with {cube.coord_system()}")

    sect.attributes['sect_info'] = sect_info
    sect.attributes['sect_info']['pnts'] = pnts
    return sect


def _add_real_lonlat(cube, lon_val, lat_val):
    geogcs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
    lon_coord = iris.coords.AuxCoord(lon_val,
                                     standard_name='longitude',
                                     units='degrees_east',
                                     coord_system=geogcs)
    lat_coord = iris.coords.AuxCoord(lat_val,
                                     standard_name='latitude',
                                     units='degrees_north',
                                     coord_system=geogcs)
    cube.add_aux_coord(lon_coord, 1)
    cube.add_aux_coord(lat_coord, 1)
