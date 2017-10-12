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
        if isinstance(out, (tuple, list, set)):
            res = []
            for iout in out:
                ires = Cube(iout,
                            dim_coords_and_dims=[(c, a_cube.coord_dims(c))
                                                 for c in a_cube.dim_coords],
                            aux_coords_and_dims=[(c, a_cube.coord_dims(c))
                                                 for c in a_cube.aux_coords],
                            units=(iout.units.__str__()
                                   .replace(' ** ', '^')
                                   .replace(' * ', ' ')))
                res.append(ires)
            res = tuple(res)
        else:
            res = Cube(out,
                       dim_coords_and_dims=[(c, a_cube.coord_dims(c))
                                            for c in a_cube.dim_coords],
                       aux_coords_and_dims=[(c, a_cube.coord_dims(c))
                                            for c in a_cube.aux_coords],
                       units=(out.units.__str__()
                              .replace(' ** ', '^')
                              .replace(' * ', ' ')))
        return res
    return wrapper


coriolis_parameter = cubehandler(calc.coriolis_parameter)
potential_temperature = cubehandler(calc.potential_temperature)
equivalent_potential_temperature = cubehandler(calc.equivalent_potential_temperature)  # noqa
vapor_pressure = cubehandler(calc.vapor_pressure)
dewpoint = cubehandler(calc.dewpoint)


def cape_and_cin(pressure, temperature, specific_humidity):
    r"""Calculate CAPE and CIN from 1-D profiles of pressure,
    temperature, and specific humidity

    Calculates parcel's profile internally, starting from the highest pressure

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
    `iris.cube.Cube`, scalar
        CAPE (J/kg)
    `iris.cube.Cube`, scalar
        CIN (J/kg)
    """
    if ([i.ndim != 1 for i in (pressure,
                               temperature,
                               specific_humidity)]).any():
        raise NotImplementedError('input cubes should be 1D')
    mixr = specific_humidity_to_mixing_ratio(specific_humidity)
    e = vapor_pressure(pressure, mixr)
    tdew = dewpoint(e)
    # p = pres.data * metunits.units(str(pressure.units))
    # t = temperature.data * metunits.units(str(temperature.units))
    # td = tdew.data * metunits.units(str(tdew.units))
    pprof = cubehandler(calc.parcel_profile)(pressure, temperature[0], tdew[0])
    cape, cin = cubehandler(calc.cape_cin)(pressure, temperature, tdew, pprof)
    return cape, cin


def specific_humidity_to_mixing_ratio(specific_humidity):
    r"""Calculate mixing ratio from specific humidity
    Parameters
    ----------
    specific_humidity : `iris.cube.Cube`
        dimensionless (kg/kg) atmospheric specific humidity
    Returns
    -------
    `iris.cube.Cube`
        Mixing ratio (kg/kg)
    """
    spechum = specific_humidity.copy()
    spechum.convert_units(1)
    mixr = spechum / ((spechum)*(-1) + 1)
    return mixr


def specific_to_relative_humidity(pressure, temperature, specific_humidity):
    r"""Calculate air relative humidity from specific humidity
    Given total `pressure` and `specific_humidity`, calculates the
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
    mixr = specific_humidity_to_mixing_ratio(specific_humidity)
    e = vapor_pressure(pressure, mixr)
    relh = e / e_s
    relh.rename('relative_humidity')
    return relh
