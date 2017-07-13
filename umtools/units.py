# -*- coding: utf-8 -*-
"""
Units
"""
import cf_units


class GrdStep(object):
    def __init__(self, key):
        assert isinstance(key, str)
        self.key = key

    def __str__(self):
        return self.key

    def to_str(self, mask='{0}.{1}km'):
        """km2p2 -> 2.2km"""
        return mask.format(self.key[2], self.key[4])

    def to_flt(self, unit='m'):
        """km2p2 -> 2200"""
        un = cf_units.Unit(self.to_str())
        return un.convert(1, cf_units.Unit(unit))


def convert_unit_str(str1, str2):
    return cf_units.Unit(str1).convert(1, cf_units.Unit(str2))
