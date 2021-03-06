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

from .grid import unrotate_lonlat_grids
from .cart import lcc_map_grid


def make_axesgrid(vrbls_lists, axsize=8, axes_pad=0.2):
    """
    TODO: deprecate?
    """
    nplots = len(vrbls_lists)
    nrows = int(np.sqrt(nplots))
    ncols = int(np.ceil(nplots / nrows))
    fig = plt.figure(figsize=(ncols*axsize, nrows*axsize))
    vrbls = vrbls_lists[0]
    # Check if colorbar is needed
    axgr_kw = dict(aspect=False, axes_pad=axes_pad)
    for icube in vrbls:
        if isinstance(icube, iris.cube.Cube):
            cbar = icube.attributes.get('colorbar')
            if cbar or isinstance(cbar, dict):
                axgr_kw.update(cbar_location='right',
                               cbar_mode='single',
                               cbar_pad=0.1,
                               cbar_size='3%')
                break
    axgr = AxesGrid(fig, 111, (nrows, ncols), **axgr_kw)
    return fig, axgr


def prepare_map(vrbls_lists, geoax=False, axsize=8):
    """
    Prepare grid of subplots for a list of lists of cubes

    Parameters
    ----------
    vrbls_lists: list of lists of cubes
        List of lists or cubelists each of which is plotted in one subplot.
        For example, it can be a list of model experiments output:
        [
          [<iris 'Cube' of air_pressure / (Pa) (model_level_number: 40; grid_latitude: 600; grid_longitude: 600)>,
           <iris 'Cube' of air_potential_temperature / (K) (model_level_number: 40; grid_latitude: 600; grid_longitude: 600)>,
           <iris 'Cube' of geopotential_height / (m) (model_level_number: 40; grid_latitude: 600; grid_longitude: 600)>],
          [<iris 'Cube' of air_pressure / (Pa) (model_level_number: 40; grid_latitude: 600; grid_longitude: 600)>,
           <iris 'Cube' of air_potential_temperature / (K) (model_level_number: 40; grid_latitude: 600; grid_longitude: 600)>,
           <iris 'Cube' of geopotential_height / (m) (model_level_number: 40; grid_latitude: 600; grid_longitude: 600)>]
        ],
        For each of the 2 lists here a subplot will be created and 3 variables will be plotted. 
    geoax: bool, optional (default False)
        Make subplots `cartopy.mpl.geoaxes.GeoAxes` (see `lcc_map_grid()`)
    axsize: float, optional (default 8)
        Subplot size in inches
    Returns
    -------
    Figure, AxesGrid, dict
        Parent figure and grid of subplots
        and a dictionary with transform keyword if geoax is True
    """
    nplots = len(vrbls_lists)
    nrows = int(np.sqrt(nplots))
    ncols = int(np.ceil(nplots / nrows))
    fig = plt.figure(figsize=(ncols*axsize, nrows*axsize))
    vrbls = vrbls_lists[0]  # TODO: allow different domains
    # Check if colorbar is needed
    axgr_kw = dict(axes_pad=0.1)
    for icube in vrbls:
        if isinstance(icube, iris.cube.Cube):
            cbar = icube.attributes.get('colorbar')
            if cbar or isinstance(cbar, dict):
                axgr_kw.update(cbar_location='right',
                               cbar_mode='single',
                               cbar_pad=0.05,
                               cbar_size='3%')
                break
    if geoax:
        lon1, lon2, lat1, lat2 = [[] for _ in range(4)]
        for icube in vrbls:
            if isinstance(icube, (list, tuple)):
                icube = icube[0]
            lons, lats = unrotate_lonlat_grids(icube)
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


def add_xaxis_below(parent_ax, xtick_array, xlab_array, shift_down):
    """
    Add x-axis parallel to the parent x-axis

    Parameters
    ----------
    parent_ax: matplotlib.axes._subplots.AxesSubplot
        Parent axes
    xtick_array: numpy.ndarray or list
        Sequence of x-tick positions
    xlab_array: numpy.ndarray or list
        Sequence of labels for x-ticks
    shift_down: int
        Number of points to shift the axis down
    Returns
    -------
    matplotlib.axes._subplots.AxesSubplot
        Added x-axis
    """
    newax = parent_ax.twiny()
    newax.set_xticks(xtick_array)
    newax.set_xticklabels(xlab_array)
    newax.spines['left'].set_visible(False)
    newax.spines['right'].set_visible(False)
    newax.set_frame_on(True)
    newax.patch.set_visible(False)
    newax.xaxis.set_ticks_position('bottom')
    newax.xaxis.set_label_position('bottom')
    newax.spines['bottom'].set_position(('outward', shift_down))
    newax.tick_params(axis='both', which='major')
    newax.grid('off')
    return newax
