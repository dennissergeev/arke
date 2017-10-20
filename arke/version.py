from os.path import join as pjoin

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
_version_major = 0
_version_minor = 4
_version_micro = 2  # use '' for first of series, number for 1 and above
_version_extra = ''
#_version_extra = ''  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro:
    _ver.append(_version_micro)
if _version_extra:
    _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

CLASSIFIERS = ["Development Status :: 3 - Alpha",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: MIT License",
               "Operating System :: POSIX :: Linux",
               "Programming Language :: Python :: 3.5",
               "Programming Language :: Python :: 3.6",
               "Topic :: Scientific/Engineering :: Atmospheric Science"]

# Description should be a one-liner:
description = "arke: a Python toolbox for UK Met Office Unified Model output"
# Long description will go up on the pypi page
# long_description = """
# 
# """

NAME = "arke"
MAINTAINER = "Denis Sergeev"
MAINTAINER_EMAIL = "dennis.sergeev@gmail.com"
DESCRIPTION = description
# LONG_DESCRIPTION = long_description
URL = "http://github.com/dennissergeev/arke"
DOWNLOAD_URL = ""
LICENSE = "MIT"
AUTHOR = "Denis Sergeev"
AUTHOR_EMAIL = "dennis.sergeev@gmail.com"
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
PACKAGES = ['arke',
            'arke.tests']
PACKAGE_DATA = {'arke': [pjoin('data', '*')]}
PYTHON_REQUIRES='>=3'
