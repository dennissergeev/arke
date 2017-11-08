# -*- coding: utf-8 -*-
"""
Meteorological constants as iris cubes

Taken from metpy.constants submodule
"""
import iris
import metpy.constants as metconst
import metpy.units as metunits

namespace = locals()

OPT_LEN = 7  # optimal-ish length of unit name

quants = {k: v for k, v in metconst.__dict__.items()
          if isinstance(v, metunits.units.Quantity)}

mag_units = set([(i.magnitude, i.units) for i in quants.values()])

nonrepeated = dict()
for k, v in quants.items():
    if len(k) <= OPT_LEN:
        nonrepeated[k] = v

counts = {k: 0 for k in nonrepeated.keys()}
for k, v in quants.items():
    key = [i for i, j in nonrepeated.items() if j == v][0]
    if (v.magnitude, v.units) in mag_units:
        counts[key] += 1

for k, v in nonrepeated.items():
    if v.units.dimensionless:
        un = 1
    else:
        un = v.units.__str__().replace(' ** ', '^').replace(' * ', ' ')
    if counts[k] == 1:
        namespace[k] = iris.cube.Cube(v.magnitude,
                                      long_name=k,
                                      units=un)
    else:
        long_name = [i for i, j in quants.items() if j == v and i != k][0]
        namespace[k] = iris.cube.Cube(v.magnitude,
                                      long_name=long_name,
                                      units=un)


namespace['Rd'].convert_units('joule / kg / kelvin')

namespace['egr_factor'] = 0.31
