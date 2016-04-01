# -*- coding: utf-8 -*-
import os
import time

import cf_units
import iris

iris.FUTURE.netcdf_no_unlimited = True

def unrotate_xy_grids(cube):
    """
    Convert rotated-pole lons and lats to unrotated ones using X and Y coordinate for a given cube.
    Unites get_xy_grids() and unrotate_pole() functions.

    Args:

        * cube - The cube with rotated coordinate system for which to generate 2D X and Y unrotated coordinates.

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


def nearest_plevel(cube, pval):
    pcoord = cube.coord('pressure')
    i = pcoord.nearest_neighbour_index(pval)
    return pcoord.points[i]


def nearest_tval(cube, dt):
    timevar = cube.coord('time')
    itime = timevar.nearest_neighbour_index(timevar.units.date2num(dt))
    return timevar.units.num2date(timevar.points[itime])

class grdstep(object):
    def __init__(self, key):
        self.key = key

    def to_str(self, mask='{0}.{1}km'):
        """km2p2 -> 2.2km"""
        return mask.format(self.key[2], self.key[4])

    def to_flt(self, unit='m'):
        """km2p2 -> 2200"""
        un = cf_units.Unit(self.to_str())
        return un.convert(1, cf_units.Unit(unit))


def stash_id_to_name(stash_id, path_to_stash=None, default_name='unknown'):
    """ Retrieve a name from the STASH table """
    try:
        import pandas as pd
    except ImportError:
        print('Unable to import pandas, default name returned instead')
        return default_name

    url = 'http://reference.metoffice.gov.uk/um/stash?_format=csv&_view=with_metadata'
    if path_to_stash is None:
        path_to_stash = os.path.join(os.curdir,'stash.csv')
        try:
            df = pd.read_csv(path_to_stash)
        except:
            print('File stash.csv not found')
            print('Trying to download it ...')
            try:
                import urllib
                f = urllib.request.URLopener()
                f.retrieve(url, path_to_stash)
            except:
                print('Download failed')
                print('Default name returned instead')

                return default_name

    df = pd.read_csv(path_to_stash)

    stash_label = df['rdfs:label'][df['@notation']==stash_id]
    if len(stash_label) > 0:
        return stash_label.values[0]
    else:
        print('Match not found, default name returned instead')
        return default_name


def convert_unit_str(str1, str2):
    return cf_units.Unit(str1).convert(1, cf_units.Unit(str2))


# tic toc functions ala Matlab
# Source: https://gist.github.com/tylerhartley/5174230
def tic(tag=None):
    '''Start timer function.
    tag = used to link a tic to a later toc. Can be any dictionary-able key.
    '''
    global TIC_TIME
    if tag is None:
        tag = 'default'

    try:
        TIC_TIME[tag] = time.time()
    except NameError:
        TIC_TIME = {tag: time.time()}


def toc(tag=None, save=False, fmt=False):
    '''Timer ending function.
    tag - used to link a toc to a previous tic. Allows multipler timers, nesting timers.
    save - if True, returns float time to out (in seconds)
    fmt - if True, formats time in H:M:S, if False just seconds.
    '''
    global TOC_TIME
    template = 'Elapsed time is:'
    if tag is None:
        tag = 'default'
    else:
        template = '%s - '%tag + template

    try:
        TOC_TIME[tag] = time.time()
    except NameError:
        TOC_TIME = {tag: time.time()}

    if TIC_TIME:
        d = (TOC_TIME[tag]-TIC_TIME[tag])

        if fmt:
            print(template + ' %s'%time.strftime('%H:%M:%S', time.gmtime(d)))
        else:
            print(template + ' %f seconds'%(d))

        if save: return d

    else:
        print("no tic() start time available. Check global var settings")
