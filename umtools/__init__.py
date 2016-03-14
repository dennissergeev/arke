from . import irismode
from .version import __version__

__all__ = []

try:
    from . import xraymode
    __all__.append('xraymode')
except ImportError:
    pass
