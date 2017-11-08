# -*- coding: utf-8 -*-
import os

import iris

iris.FUTURE.netcdf_no_unlimited = True


def stash_id_to_name(stash_id, path_to_stash=None, default_name='unknown'):
    """ Retrieve a name from the STASH table """
    try:
        import pandas as pd
    except ImportError:
        print('Unable to import pandas, default name returned instead')
        return default_name

    url = ('http://reference.metoffice.gov.uk'
           '/um/stash?_format=csv&_view=with_metadata')
    if path_to_stash is None:
        path_to_stash = os.path.join(os.curdir, 'stash.csv')
        try:
            df = pd.read_csv(path_to_stash)
        except (FileNotFoundError, pd.errors.ParserError):
            print('File stash.csv not found')
            print('Trying to download it ...')
            try:
                import urllib
                f = urllib.request.URLopener()
                f.retrieve(url, path_to_stash)
            except urllib.error.HTTPError:
                print('Download failed')
                print('Default name returned instead')
                return default_name

    df = pd.read_csv(path_to_stash)

    stash_label = df['rdfs:label'][df['@notation'] == stash_id]
    if len(stash_label) > 0:
        return stash_label.values[0]
    else:
        print('Match not found, default name returned instead')
        return default_name


def replace_unknown_names(dataset, default_name='unknown'):
    """
    Replace (in place) missing names within an `iris.cube.CubeList`
    using STASH code
    """
    for ivar in dataset:
        if default_name in ivar.name().lower():
            try:
                stash_id = ivar.attributes['STASH'].__str__()
                ivar.rename(stash_id_to_name(stash_id))
            except AttributeError:
                print('Unable to rename, STASH attribute is missing')


# def rename_cubes_from_stashmaster(cube, stashmaster_file,
#                                   parser=configparser.ConfigParser(),
#                                   add_help=True):
#     if cube.name().lower() == 'unknown':
#         stash = cube.attributes['STASH']
#         parser.read(stashmaster_file)
#         try:
#             section = parser['stashmaster:code({})'.format(stash.item)]
#             new_name = '_'.join(section['description'].split())
#             cube.rename(new_name)
#             if add_help:
#                 try:
#                     cube.attributes['help'] = section['help'].\
#                                               replace('=', '').\
#                                               replace('\n', '')
#                 except KeyError:
#                     pass
#         except KeyError:
#             pass
