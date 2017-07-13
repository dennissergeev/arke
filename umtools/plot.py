# -*- coding: utf-8 -*-
"""
Plot hourly MetUM output on the same figure
"""
# Standard packages
import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np

from .numerics import unrotate_xy_grids
from .cart import lcc_map_grid


def prepare_map(vrbls_lists, geoax=False):
    """
    TODO: docstring
    """
    nplots = len(vrbls_lists)
    nrows = int(np.sqrt(nplots))
    ncols = int(np.ceil(nplots / nrows))
    fig = plt.figure(figsize=(ncols*8, nrows*8))
    vrbls = vrbls_lists[0]  # TODO: allow different domains
    # Check if colorbar is needed
    axgr_kw = {}
    for icube in vrbls:
        if isinstance(icube, iris.cube.Cube):
            cbar = icube.attributes.get('colorbar')
            if cbar or isinstance(cbar, dict):
                axgr_kw.update(axes_pad=0.1,
                               cbar_location='right',
                               cbar_mode='single',
                               cbar_pad=0.05,
                               cbar_size='3%')
                break
    if geoax:
        lon1, lon2, lat1, lat2 = [[] for _ in range(4)]
        for icube in vrbls:
            if isinstance(icube, list):
                icube = icube[0]
            lons, lats = unrotate_xy_grids(icube)
            if hasattr(icube.data, 'mask'):
                _mask = icube.data.mask
                lons = np.ma.masked_where(_mask, lons)
                lats = np.ma.masked_where(_mask, lats)
                lon1.append(np.min(lons))
                lon2.append(np.max(lons))
            else:
                lon1.append(np.mean(lons[:,  0]))
                lon2.append(np.mean(lons[:, -1]))
            lat1.append(np.min(lats))
            lat2.append(np.max(lats))
        lon1 = np.min(lon1)
        lon2 = np.max(lon2)
        lat1 = np.min(lat1)
        lat2 = np.max(lat2)
        # xtick = best_ticks[np.argmin(abs(np.array(best_ticks)
        #                                  - (lon2-lon1) * 0.1))]
        # ytick = best_ticks[np.argmin(abs(np.array(best_ticks)
        #                                  - (lat2-lat1) * 0.5))]
        ticks = None  # [xtick, ytick]
        clon = 0.5 * (lon1 + lon2)
        clat = 0.5 * (lat1 + lat2)
        extent = [lon1, lon2, lat1, lat2]
        # coast = dict(scale='50m', edgecolor='#AAAAAA', facecolor='#FFFFFF')
        coast = dict(scale='50m', edgecolor='0.75', alpha=0.5, facecolor='0.5')
        lcc_kw = dict(clon=clon, clat=clat, coast=coast,
                      extent=extent, ticks=ticks)
        axgr = lcc_map_grid(fig, (nrows, ncols), **lcc_kw, **axgr_kw)
        # cax = axgr.cbar_axes[0]
        mapkey = dict(transform=ccrs.PlateCarree())
    else:
        axgr = AxesGrid(fig, 111, (nrows, ncols), **axgr_kw)
        # cax = axgr.cbar_axes[0]
        # ax = fig.add_subplot(111)
        mapkey = {}

    return fig, axgr, mapkey
