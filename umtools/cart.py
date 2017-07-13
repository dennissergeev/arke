# -*- coding: utf-8 -*-
"""
Collection of functions for cartographic plotting.
Based on cartopy.
"""
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
from copy import copy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import shapely.geometry as sgeom

best_ticks = [1, 2, 5, 10, 15, 20, 30, 50]


def find_side(ls, side):
    """
    Given a shapely LineString which is assumed to be rectangular, return the
    line corresponding to a given side of the rectangle.

    """
    minx, miny, maxx, maxy = ls.bounds
    points = {'left': [(minx, miny), (minx, maxy)],
              'right': [(maxx, miny), (maxx, maxy)],
              'bottom': [(minx, miny), (maxx, miny)],
              'top': [(minx, maxy), (maxx, maxy)], }
    return sgeom.LineString(points[side])


def lambert_xticks(ax, ticks):
    """Draw ticks on the bottom x-axis of a Lambert Conformal projection."""
    def _te(xy):
        return xy[0]

    def _lc(t, n, b):
        return np.vstack((np.zeros(n) + t, np.linspace(b[2], b[3], n))).T
    xticks, xticklabels = _lambert_ticks(ax, ticks, 'bottom', _lc, _te)
    ax.xaxis.tick_bottom()
    ax.set_xticks(xticks)
    ax.set_xticklabels([ax.xaxis.get_major_formatter()(xtick)
                        for xtick in xticklabels])


def lambert_yticks(ax, ticks):
    """Draw ticks on the left y-axis of a Lamber Conformal projection."""
    def _te(xy):
        return xy[1]

    def _lc(t, n, b):
        return np.vstack((np.linspace(b[0], b[1], n), np.zeros(n) + t)).T
    yticks, yticklabels = _lambert_ticks(ax, ticks, 'left', _lc, _te)
    ax.yaxis.tick_left()
    ax.set_yticks(yticks)
    ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick)
                        for ytick in yticklabels])


def _lambert_ticks(ax, ticks, tick_location, line_constructor, tick_extractor):
    """
    Get the tick locations and labels for
    an axis of a Lambert Conformal projection.
    """
    outline_patch = sgeom.LineString(ax.outline_patch.get_path().
                                     vertices.tolist())
    axis = find_side(outline_patch, tick_location)
    n_steps = 30
    extent = ax.get_extent(ccrs.PlateCarree())
    _ticks = []
    for t in ticks:
        xy = line_constructor(t, n_steps, extent)
        proj_xyz = ax.projection.transform_points(ccrs.Geodetic(),
                                                  xy[:, 0], xy[:, 1])
        xyt = proj_xyz[..., :2]
        ls = sgeom.LineString(xyt.tolist())
        locs = axis.intersection(ls)
        if not locs:
            tick = [None]
        else:
            tick = tick_extractor(locs.xy)
        _ticks.append(tick[0])
    # Remove ticks that aren't visible:
    ticklabels = copy(ticks)
    while True:
        try:
            index = _ticks.index(None)
        except ValueError:
            break
        _ticks.pop(index)
        ticklabels.pop(index)
    return _ticks, ticklabels


def label_map(ax, toponyms={}, transform=None):
    if transform is None:
        transform = ax.transAxes
    for i in toponyms:
        txt = ax.text(i['lon'], i['lat'], i['name'],
                      transform=transform,
                      color='k', fontsize=20, fontweight='bold',
                      ha='center', va='center')
        txt.set_zorder(200)


def pc_map(fig, subplot_grd=111,
           projection=ccrs.PlateCarree(), coast=None, extent=None):
    ax = fig.add_subplot(subplot_grd, projection=projection)

    gl = ax.gridlines(draw_labels=True)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    if isinstance(extent, list):
        ax.set_extent(extent, crs=ccrs.PlateCarree())
    add_coastline(ax, coast)
    return ax


def get_xy_ticks(ticks):
    # Define gridline locations
    try:
        if len(ticks) == 2:
            try:
                _x, _y = ticks[0], ticks[1]
            except KeyError:
                _x, _y = ticks['x'], ticks['y']

            if (isinstance(_x, (int, float))
               and isinstance(_y, (int, float))):
                # assume ticks is [xstep, ystep] sequence
                xticks = list(np.arange(-180, 181, _x))
                yticks = list(np.arange(-90, 91, _y))
            elif (isinstance(_x, (tuple, list, np.ndarray))
                  and isinstance(_y, (tuple, list, np.ndarray))):
                # assume ticks is [xticks, yticks] sequence
                xticks, yticks = list(_x), list(_y)
    except TypeError:
        # fall back to default arrays
        xticks = list(np.linspace(-180, 180, 37))
        yticks = list(np.linspace(-90, 90, 19))
    return xticks, yticks


def lcc_map(fig, subplot_grd=111, clon=None, clat=None, extent=None,
            coast=None, ticks=None):
    # Create a Lambert Conformal projection:
    proj = ccrs.LambertConformal(central_longitude=clon, central_latitude=clat)

    # Draw a set of axes with coastlines:
    ax = fig.add_subplot(subplot_grd, projection=proj)
    if isinstance(extent, list):
        ax.set_extent(extent, crs=ccrs.PlateCarree())

    add_coastline(ax, coast)

    if ticks:
        xticks, yticks = get_xy_ticks(ticks)
        # *must* call draw in order to get the axis boundary used to add ticks
        fig.canvas.draw()
        # Draw the lines using cartopy's built-in gridliner
        ax.gridlines(xlocs=xticks, ylocs=yticks)
        # Label the end-points of the gridlines using the custom tick makers:
        ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
        ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
        lambert_xticks(ax, xticks)
        lambert_yticks(ax, yticks)

    return ax


# soon to be included in cartopy
class GeoAxesGrid(AxesGrid):
    def __init__(self, fig, rect, nrows_ncols, projection, **axesgrid_kw):
        axes_class = (GeoAxes, dict(map_projection=projection))
        axesgrid_kw['label_mode'] = ''
        super(GeoAxesGrid, self).__init__(fig, rect, nrows_ncols,
                                          axes_class=axes_class, **axesgrid_kw)


def lcc_map_grid(fig, nrows_ncols, clon=None, clat=None, extent=None,
                 coast=None, ticks=None, **axesgrid_kw):
    proj = ccrs.LambertConformal(central_longitude=clon, central_latitude=clat)
    axgr = GeoAxesGrid(fig, 111, nrows_ncols, projection=proj, **axesgrid_kw)
    if ticks is not None:
        xticks, yticks = get_xy_ticks(ticks)

    for ax in axgr:
        if isinstance(extent, list):
            ax.set_extent(extent, crs=ccrs.PlateCarree())

        add_coastline(ax, coast)

        if ticks is not None:
            # *must* call draw in order to get the axis boundary to add ticks
            fig.canvas.draw()
            # Draw the lines using cartopy's built-in gridliner
            ax.gridlines(xlocs=xticks, ylocs=yticks)
            # Label the end-points of the gridlines using the tick makers:
            ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
            ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
            lambert_xticks(ax, xticks)
            lambert_yticks(ax, yticks)

    return axgr


def add_coastline(ax, coast):
    if isinstance(coast, str):
        feature = cartopy.feature.NaturalEarthFeature(name='coastline',
                                                      category='physical',
                                                      scale=coast,
                                                      edgecolor='#AAAAAA',
                                                      facecolor='#AAAAAA')
        ax.add_feature(feature)
    elif isinstance(coast, dict):
        feature = cartopy.feature.NaturalEarthFeature(name='coastline',
                                                      category='physical',
                                                      **coast)
        ax.add_feature(feature)
