# -*- coding: utf-8 -*-
"""
Collection of functions for cartographic plotting.
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


def _find_side(ls, side):
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


def _lambert_xticks(ax, ticks):
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


def _lambert_yticks(ax, ticks):
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
    axis = _find_side(outline_patch, tick_location)
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


def label_map(ax, toponyms, transform=None, **text_kw):
    """
    Put labels on a map

    Parameters
    ----------
    ax: cartopy.mpl.geoaxes.GeoAxesSubplot
        axes to put names on
    toponyms: list
        list of dictionaries containing `lon`, `lat`, `name` keys
        defining the longitude and latitude of each `name` toponym
    transform: matplotlib.transforms.BboxTransformTo, optional
        axes transform; set to ax.transAxes by default
    """
    if transform is None:
        transform = ax.transAxes
    for i in toponyms:
        txt = ax.text(i['lon'], i['lat'], i['name'],
                      transform=transform,
                      **text_kw)
        txt.set_zorder(20)


def pc_map(fig, subplot_grd=111,
           projection=ccrs.PlateCarree(), coast=None, extent=None):
    """
    Create axes with the Plate Carree projection in a given figure

    Parameters
    ----------
    fig: matplotlib.figure.Figure
         matplotlib figure
    subplot_grd: int, optional
        3-digit integer describing the position of the subplot
        default: 111
    projection: str or cartopy.crs.CRS, optional
        projection class of the axes, default: cartopy.crs.PlateCarree
    coast: str or dict, optional
        parameters to draw a coastline, see `add_coastline()` for details
    extent: sequence, optional
        extent (x0, x1, y0, y1) of the map in the given coordinate projection
    Returns
    -------
    cartopy.mpl.geoaxes.GeoAxesSubplot
        axes with the Plate Caree projection
    """
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
    """Define gridline locations"""
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
    """
    Create axes the Lambert Conformal Conic (LCC) in a given figure

    Parameters
    ----------
    fig: matplotlib.figure.Figure
         matplotlib figure
    subplot_grd: int, optional
        3-digit integer describing the position of the subplot
        default: 111
    clon: float
        central longitude of LCC projection
    clat: float
        central latitude of LCC projection
    coast: str or dict, optional
        parameters to draw a coastline, see `add_coastline()` for details
    extent: sequence, optional
        extent (x0, x1, y0, y1) of the map in the given coordinate projection
    ticks: sequence, optional
        see `get_xy_ticks()` for details
    Returns
    -------
    cartopy.mpl.geoaxes.GeoAxesSubplot
        axes with the LCC projection
    """
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
        _lambert_xticks(ax, xticks)
        _lambert_yticks(ax, yticks)

    return ax


class GeoAxesGrid(AxesGrid):
    """
    A subclass of :class:`mpl_toolkits.axes_grid1.AxesGrid` representing
    a grid of maps with the same projection :class:`~cartopy.crs.Projection`.
    .. note::
       * `axes_class` is defined automatically
       * The :class:`AxesGrid` built-in labelling is always switched off,
         and instead a standard procedure of creating
         grid lines and labels should be used.
    """
    def __init__(self, fig, rect, nrows_ncols, projection, **axesgrid_kw):
        """
        Build a :class:`GeoAxesGrid` instance with a grid nrows*ncols
        :class:`GeoAxes` with a projection :class:`~cartopy.crs.Projection`
        in :class:`~matplotlib.figure.Figure` *fig* with
        *rect=[left, bottom, width, height]* (in
        :class:`~matplotlib.figure.Figure` coordinates) or
        the subplot position code (e.g., "121").
        Kwargs:
          Keyword           Default   Description
          ================  ========  =========================================
          direction         "row"     [ "row" | "column" ]
          axes_pad          0.02      float| pad between axes given in inches
                                      or tuple-like of floats,
                                      (horizontal padding, vertical padding)
          add_all           True      [ True | False ]
          share_all         False     [ True | False ]
          aspect            True      [ True | False ]
          cbar_mode         None      [ "each" | "single" | "edge" ]
          cbar_location     "right"   [ "left" | "right" | "bottom" | "top" ]
          cbar_pad          None
          cbar_size         "5%"
          cbar_set_cax      True      [ True | False ]
          ================  ========  =========================================
        *cbar_set_cax* : if True, each axes in the grid has a cax
          attribute that is bind to associated cbar_axes.
        """
        axesgrid_kw['axes_class'] = (GeoAxes, dict(map_projection=projection))
        axesgrid_kw['label_mode'] = ''  # note the empty label_mode
        super(GeoAxesGrid, self).__init__(fig, rect, nrows_ncols,
                                          **axesgrid_kw)


def lcc_map_grid(fig, nrows_ncols, clon, clat, extent=None,
                 coast=None, ticks=None, **axesgrid_kw):
    """
    Build an `AxesGrid` instance with a grid `nrows_ncols` with
    the Lambert Conformal Conic (LCC) projection and `**axesgrid_kw` parameters

    Parameters
    ----------
    fig: matplotlib.figure.Figure
        parent figure
    nrows_ncols: tuple of int
        N rows and N cols
    clon: float
        central longitude of LCC projection
    clat: float
        central latitude of LCC projection
    coast: str or dict, optional
        parameters to draw a coastline, see `add_coastline()` for details
    extent: sequence, optional
        extent (x0, x1, y0, y1) of the map in the given coordinate projection
    ticks: sequence, optional
        see `get_xy_ticks()` for details
    Returns
    -------
    GeoAxesGrid
    """
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
            _lambert_xticks(ax, xticks)
            _lambert_yticks(ax, yticks)

    return axgr


def merc_map_grid(fig, nrows_ncols, extent=None,
                  coast=None, ticks=None, **axesgrid_kw):
    """
    Build an `AxesGrid` instance with a grid `nrows_ncols` with
    the Mercator projection and `**axesgrid_kw` parameters

    Parameters
    ----------
    fig: matplotlib.figure.Figure
        parent figure
    nrows_ncols: tuple of int
        N rows and N cols
    coast: str or dict, optional
        parameters to draw a coastline, see `add_coastline()` for details
    extent: sequence, optional
        extent (x0, x1, y0, y1) of the map in the given coordinate projection
    ticks: sequence, optional
        see `get_xy_ticks()` for details
    axesgrid_kw: dict, optional
        AxesGrid class keywords
    Returns
    -------
    GeoAxesGrid
    """
    proj = ccrs.Mercator()
    axgr = GeoAxesGrid(fig, 111, nrows_ncols, projection=proj, **axesgrid_kw)
    if ticks is not None:
        xticks, yticks = get_xy_ticks(ticks)

    for ax in axgr:
        if isinstance(extent, list):
            ax.set_extent(extent, crs=ccrs.PlateCarree())

        add_coastline(ax, coast)

        if ticks is not None:
            ax.gridlines(xlocs=xticks, ylocs=yticks)
            ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
            ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)

    return axgr


def add_coastline(ax, coast):
    """
    Add coast outline to a given GeoAxes

    Parameters
    ----------
    ax: cartopy.mpl.geoaxes.GeoAxesSubplot
        axes to add coastlines
    coast: str or dict
        If str object is given it assumed to be a named `scale`
        (resolution to use from the Natural Earth dataset,
         currently can be one of "110m", "50m", and "10m").
        If dict is given, it assumed to contain the scale, as well as
        other kwargs of `cartopy.feature.NaturalEarthFeature`.
    """
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
