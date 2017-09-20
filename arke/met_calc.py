# -*- coding: utf-8 -*-
"""
Meteorological calculations

Taken from metpy.calc submodule and decorated to handle iris cubes
"""
from functools import wraps

from iris.cube import Cube
from metpy import calc
import metpy.units as metunits


def cubehandler(f):
    @wraps(f)
    def wrapper(*args, **kwds):
        nargs = []
        for arg in args:
            if isinstance(arg, Cube):
                a_cube = arg
                q = arg.data * metunits.units(str(arg.units))
                nargs.append(q)
            else:
                nargs.append(arg)
        out = f(*nargs, **kwds)
        cube_out = Cube(out,
                        dim_coords_and_dims=[(c, a_cube.coord_dims(c))
                                             for c in a_cube.dim_coords],
                        aux_coords_and_dims=[(c, a_cube.coord_dims(c))
                                             for c in a_cube.aux_coords],
                        units=(out.units.__str__()
                               .replace(' ** ', '^')
                               .replace(' * ', ' ')))
        return cube_out
    return wrapper


coriolis_parameter = cubehandler(calc.coriolis_parameter)
potential_temperature = cubehandler(calc.potential_temperature)
equivalent_potential_temperature = cubehandler(calc.equivalent_potential_temperature)  # noqa
vapor_pressure = cubehandler(calc.vapor_pressure)
dewpoint = cubehandler(calc.dewpoint)


def specific_to_relative_humidity(pressure, temperature, specific_humidity):
    r"""Calculate air relative humdity from specific humidity.
    Given total `pressure` and `specific_humdidity`, calculates the
    relative humidity.

    Parameters
    ----------
    pressure : `iris.cube.Cube`
        total atmospheric pressure
    temperature : `iris.cube.Cube`
        atmospheric temperature
    specific_humidity : `iris.cube.Cube`
        dimensionless (kg/kg) atmospheric specific humidity
    Returns
    -------
    `iris.cube.Cube`
        Relative humidity

    Notes
    -----
    * :math:`RH` is relative humidity
    * :math:`e` is vapor pressure calculated from the specific humidity
    * :math:`e_s` is the saturation vapor pressure
    """
    e_s = cubehandler(calc.saturation_vapor_pressure)(temperature)
    mixr = specific_humidity / ((specific_humidity)*(-1) + 1)
    e = vapor_pressure(pressure, mixr)
    relh = e / e_s
    relh.rename('relative_humidity')
    return relh
