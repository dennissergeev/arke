# -*- coding: utf-8 -*-
"""
Exceptions specific to the arke package.
"""
from iris.exceptions import ConcatenateError
from metpy.units import UndefinedUnitError


class ArkeError(Exception):
    """Base class for errors in arke package."""
    pass


class NotYetImplementedError(ArkeError):
    """
    Raised by missing functionality.

    Different meaning to NotImplementedError, which is for abstract methods.

    """
    pass


class CoordinateRankError(ArkeError):
    """Raised when coordinate has unexpected dimensionality"""
    pass


class ArgumentError(ArkeError):
    """Raised when argument type is not recognized"""
    pass


class CubeConcatenateError(ConcatenateError):
    """
    Raised when concatenate is expected to produce a single cube,
    but fails to do so.
    """
    pass


class CubeNotFoundError(ArkeError):
    """Raised when a search yields no cubes."""
    pass


class BadCoordinateError(ArkeError):
    """Raised when no suitable coordinate is found."""
    pass


class UndefinedMetUnitError(UndefinedUnitError):
    """Raised when units are undefined."""
    pass
