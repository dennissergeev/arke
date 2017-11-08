# -*- coding: utf-8 -*-
"""
Meteorological calculations

Taken from metpy.calc submodule and decorated to handle iris cubes
"""
from functools import wraps

from iris.cube import Cube
import cf_units
from metpy import calc
import metpy.units as metunits
import numpy as np


def cubehandler(f):
    """
    Wrap a function from metpy.calc submodule to pass iris cubes as arguments
    and get cubes as output.

    Copies coordinates from the first cube among the arguments.
    Only copies `dim_coords` and `aux_coords`.

    May not work with some units.
    """
    @wraps(f)
    def wrapper(*args, **kwds):
        nargs = []
        a_cube = None
        for arg in args:
            if isinstance(arg, Cube):
                if arg.ndim > 0:
                    a_cube = arg
                elif a_cube is None:
                    a_cube = arg
                for ut_format in set(cf_units.UT_FORMATS):
                    try:
                        un = metunits.units(arg.units.format(ut_format))
                    except:
                        pass
                if np.ma.is_masked(arg.data):
                    q = metunits.masked_array(arg.data, data_units=un)
                else:
                    q = arg.data * un
                nargs.append(q)
            else:
                nargs.append(arg)
        out = f(*nargs, **kwds)
        if isinstance(out, (tuple, list, set)):
            res = []
            for iout in out:
                ires = Cube(iout.magnitude,
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
            res = Cube(out.magnitude,
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
frontogenesis = cubehandler(calc.frontogenesis)


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
    if any([i.ndim != 1 for i in (pressure,
                                  temperature,
                                  specific_humidity)]):
        raise NotImplementedError('input cubes should be 1D')
    if pressure.data[0] < pressure.data[-1]:
        # reverse arrays
        zdir = slice(None, None, -1)
    else:
        zdir = slice(None)
    mixr = specific_humidity_to_mixing_ratio(specific_humidity)
    e = vapor_pressure(pressure, mixr)
    tdew = dewpoint(e)
    pres = pressure[zdir]
    temp = temperature[zdir]
    tdew = tdew[zdir]
    pprof = cubehandler(calc.parcel_profile)(pres, temp[0], tdew[0])
    cape, cin = cubehandler(calc.cape_cin)(pres, temp, tdew, pprof)
    cape.rename('atmosphere_convective_available_potential_energy')
    cin.rename('atmosphere_convective_inhibition')
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
