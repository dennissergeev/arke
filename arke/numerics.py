# -*- coding: utf-8 -*-
"""
Calculate derivatives and other useful diagnostics of atmospheric flow
"""
import numpy as np
from cached_property import cached_property
import iris
from iris.analysis.calculus import differentiate
from iris.coords import AuxCoord
from iris.cube import Cube, CubeList
import warnings

from . import met_calc as mcalc
from . import met_const as mconst
from . import grid
from . import coords

phys_coord = ('height', 'level_height',
              'atmosphere_hybrid_height_coordinate',
              'pressure')

is_physical = lambda z: z.name() in phys_coord  # NOQA
notdimless_and_1d = lambda z: not z.units.is_dimensionless() and z.ndim==1  # NOQA

iscube = lambda x: isinstance(x, Cube)  # NOQA
iscube_and_not_scalar = lambda x: (isinstance(x, Cube)  # NOQA
                                   and x.shape != ())  # NOQA


def replace_dimcoord(cube, src_cube, axes='xy', return_copy=True):
    """
    Replace `dim_coords` of one cube with `dim_coords` of another

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube that needs new dim_coords
    src_cube: iris.cube.Cube
        Cube-donor of dim_coords
    axes: str
        Coordinate axes to replace dim_coords
    return_copy: bool (default True)
        if True, replace dim_coords of a copy of the `cube`
    Returns
    -------
    (if return_copy is True) iris.cube.Cube
        copy of the given cube with new dim_coords
    """
    if return_copy:
        cp_cube = cube.copy()
    for axis in axes:
        oldcoord = cube.coord(axis=axis)
        newcoord = src_cube.coord(axis=axis)
        ndim = cube.coord_dims(oldcoord)[0]
        if return_copy:
            cp_cube.remove_coord(oldcoord)
            cp_cube.add_dim_coord(newcoord, ndim)
        else:
            cube.remove_coord(oldcoord)
            cube.add_dim_coord(newcoord, ndim)
    if return_copy:
        return cp_cube


def prepare_cube_zcoord(cube, rm_z_bounds=True, rm_z_varname=True):
    """
    Pick physically meaningful coordinate instead of level indices

    Parameters
    ----------
    cube: iris.cube.Cube
        Cube that has a z-axis coordinate
    rm_z_bounds: bool (default True)
        Remove `bounds` attribute of a selected vertical coordinate
    rm_z_varname: bool (default True)
        Remove `var_name` attribute of a selected vertical coordinate
    Returns
    -------
    iris.cube.Cube
        Copy of the original cube with replaced z-axis `dim_coord`
    """
    res = cube.copy()
    z_ax_coords = res.coords(axis='z')
    suitable_z = list(filter(notdimless_and_1d, z_ax_coords))
    if len(suitable_z) > 0:
        suitable_z = list(filter(is_physical, suitable_z))
        if len(suitable_z) == 0:
            # TODO: what if suitable_z > 1 ?
            msg = "Warning: the name of '{name}' is not among: {phys}"
            print(msg.format(name=suitable_z.name(),
                             phys=', '.join(phys_coord)))
        zcoord = suitable_z[0]
    else:
        raise ValueError('No suitable coords among: {z}'.format(z=z_ax_coords))

    if rm_z_bounds:
        zcoord.bounds = None
    if rm_z_varname:
        zcoord.var_name = None

    if zcoord not in res.dim_coords:
        for z in z_ax_coords:
            if z in res.dim_coords:
                zdim = res.coord_dims(z)[0]
                res.remove_coord(z)
        res.remove_coord(zcoord)
        zcoord = iris.coords.DimCoord.from_coord(zcoord)
        res.add_dim_coord(zcoord, zdim)
    return res


def replace_lonlat_dimcoord_with_cart(cube, dx=1, dy=None, dxdy_units='m'):
    """
    Replace latitude- and longitude-like `dim_coords`
    with points equally-spaced by `dx` grid step. Used to prepare output from
    limited-area MetUM run.

    Parameters
    ----------
    cube: iris.cube.Cube
        Target cube
    dx: float, optional (default 1)
        grid spacing of the `cube`'s x-axis
        should be in the same units as dy
    dy: float, optional (default =dx)
        grid spacing in along y-axis
    dxdy_units: str (default 'm')
        Units of grid spacings
    Returns
    -------
    iris.cube.Cube
        Copy of the original cube with cartesian x,y-coordinates
    """
    res = cube.copy()
    if dy is None:
        dy = dx
    # TODO: dx = np.diff(glats.points).mean() * 111.3
    for axis, standard_name, step in zip(('x',         'y'),
                                         ('projection_x_coordinate', 'projection_y_coordinate'),  # NOQA
                                         (dx,           dy)):
        icoord = res.coord(axis=axis)
        coord_dim = res.coord_dims(icoord)[0]
        iris.util.demote_dim_coord_to_aux_coord(res, icoord.name())

        eq_spaced_points = np.array(range(icoord.shape[0]))*step
        long_name = 'distance_along_{0}_axis'.format(axis)
        new_coord = iris.coords.DimCoord(points=eq_spaced_points,
                                         standard_name=standard_name,
                                         long_name=long_name,
                                         units=dxdy_units)

        res.add_dim_coord(new_coord, coord_dim)

    return res


def prepare_cube_on_model_levels(cube, lonlat2cart_kw={}, prep_zcoord_kw={},
                                 rm_surf_alt=True,
                                 rm_sigma=True,
                                 rm_aux_factories=True):
    """
    Prepare a cube for `AtmosFlow`:
        * replace longitude and latitude coords with cartesian coordinates
          (see `replace_lonlat_dimcoord_with_cart` for details)
        * replace vertical dim_coord with a physically-meaningful one instead
          of model levels (see `prepare_cube_zcoord` for details)

    Parameters
    ----------
    cube: iris.cube.Cube
        Target cube
    lonlat2cart_kw: dict, optional
        kwargs for `replace_lonlat_dimcoord_with_cart()`
    prep_zcoord_kw: dict, optional
        kwargs for `prepare_cube_zcoord()`
    rm_surf_alt: bool (default True)
        Remove `surface_altitude` coordinate
    rm_sigma: bool (default True)
        Remove `sigma` coordinate
    rm_aux_factories: bool (default True)
        Remove `aux_factories` of the cube
    Returns
    -------
    iris.cube.Cube
        Copy of the cube ready for building `AtmosFlow`
    """
    res = cube.copy()
    if isinstance(lonlat2cart_kw, dict):
        res = replace_lonlat_dimcoord_with_cart(res, **lonlat2cart_kw)
    if isinstance(prep_zcoord_kw, dict):
        res = prepare_cube_zcoord(res, **prep_zcoord_kw)
    # Get rid of unnecessary coordinates that hinder cube.coord() method
    if rm_surf_alt and len(res.coords('surface_altitude')) > 0:
        res.remove_coord('surface_altitude')
    if rm_sigma and len(res.coords('sigma')) > 0:
        res.remove_coord('sigma')
    if rm_aux_factories:
        # aka remove DerivedCoords
        [res.remove_aux_factory(i) for i in res.aux_factories]
    return res


def clean_pressure_coord(cube):
    """Convert pressure coordinate units to pascals and round the values"""
    try:
        pcoord = cube.coord('pressure')
        units2pa = pcoord.units.convert(1, 'Pa')
        new_pcoord = iris.coords.DimCoord((pcoord.points * units2pa).round(),
                                          long_name='pressure', units='Pa')
        dim = cube.coord_dims('pressure')[0]
        cube.remove_coord('pressure')
        cube.add_dim_coord(new_pcoord, dim)
    except iris.exceptions.CoordinateNotFoundError:
        pass


def check_coords(cubes):
    """Check the cubes coordinates for consistency"""
    # get the names of all coords binned into useful comparison groups
    coord_comparison = iris.analysis.coord_comparison(*cubes)

    bad_coords = coord_comparison['ungroupable_and_dimensioned']
    if bad_coords:
        raise ValueError("Coordinates found in one cube that describe "
                         "a data dimension which weren't in the other "
                         "cube ({}), try removing this coordinate.".format(
                             ', '.join(group.name() for group in bad_coords)))

    bad_coords = coord_comparison['resamplable']
    if bad_coords:
        raise ValueError('Some coordinates are different ({}), consider '
                         'resampling.'.format(
                             ', '.join(group.name() for group in bad_coords)))


def cube_deriv(cube, coord):
    """
    Wrapper to `iris.analysis.calculus.differentiate` to differentiate a cube
    and regrid/interpolate the result from mid-points
    to the points of the input cube
    """
    res = differentiate(cube, coord)
    kw = dict(dim_coords=True)
    if (res.coord(axis='x', **kw) != cube.coord(axis='x', **kw) or
       res.coord(axis='y', **kw) != cube.coord(axis='y', **kw)):
        res = res.regridded(cube)
    elif res.coord(axis='z', **kw) != cube.coord(axis='z', **kw):
        cube_z_points = [(cube.coord(axis='z').name(),
                          cube.coord(axis='z').points)]
        res = res.interpolate(cube_z_points, iris.analysis.Linear())
    return res


class AtmosFlow:
    """
    Atmospheric Flow

    Used to calculate meteorological parameters from the given cubes.
    Derived quantities are stored as cached properties to save computational
    time.

    Calculating quantites that involve on horizontal derivatives are only
    true if used in cartesian coordinates, e.g. on the output from a LAM model
    with constant grid spacing. Use `prepare_cube_on_model_levels` to prepare
    cubes for `AtmosFlow` on a cartesian grid.

    Attributes
    ----------
    cubes: iris.cube.CubeList
        list of cubes representing meteorological parameters
    main_cubes: iris.cube.CubeList
        list of non-scalar cubes
    wind_cmpnt: iris.cube.CubeList
        list of u,v,w-wind components
    {x,y,z}coord: iris.coord.Coord
        Coordinates in the respective dimensions
    pres: iris.cube.Cube
        Pressure created from a coordinate, if possible
    lats: iris.cube.Cube
        latitudes
    fcor: iris.cube.Cube
        Coriolis parameter (taken at 45N by default)
    d{x,y}: iris.cube.Cube
        Grid spacing (if cartesian=True)
    """
    def __init__(self, cartesian=True, **kw_vars):
        """
        Parameters
        ----------
        cartesian: bool (default True)
            Cartesian coord system flag
        **kw_vars: dict of iris cubes
            meteorological parameters

        Examples
        --------
        Initialise an `AtmosFlow` object with 3 wind components
        >>> AF = AtmosFlow(u=u_cart, v=v_cart, w=w_cart)
        and calculate relative vorticity:
        >>> rv = AF.rel_vort
        """
        self.__dict__.update(kw_vars)
        self.cartesian = cartesian

        self.cubes = CubeList(filter(iscube, self.__dict__.values()))
        self.main_cubes = CubeList(filter(iscube_and_not_scalar,
                                          self.__dict__.values()))
        self.wind_cmpnt = CubeList(filter(None,
                                          [getattr(self, 'u', None),
                                           getattr(self, 'v', None),
                                           getattr(self, 'w', None)]))
        thecube = self.main_cubes[0]

        check_coords(self.main_cubes)

        # Get the dim_coord, or None if none exist, for the xyz dimensions
        self.xcoord = thecube.coord(axis='X', dim_coords=True)
        self.ycoord = thecube.coord(axis='Y', dim_coords=True)
        self.zcoord = thecube.coord(axis='Z')
        if self.zcoord.units.is_convertible('Pa'):
            # Check if the vertical coordinate is pressure
            self.zmode = 'pcoord'
            for cube in self.main_cubes:
                if self.zcoord in cube.dim_coords:
                    clean_pressure_coord(cube)
            self.pres = coords.pres_coord_to_cube(thecube)
            self.cubes.append(self.pres)

        if not hasattr(self, 'lats'):
            try:
                _, lats = grid.unrotate_lonlat_grids(thecube)
            except (ValueError, AttributeError):
                lats = np.array([45.])
            self.lats = Cube(lats,
                             units='degrees',
                             standard_name='latitude')
        self.fcor = mcalc.coriolis_parameter(self.lats)
        self.fcor.convert_units('s-1')

        if self.cartesian:
            for ax, rot_name in zip(('x',  'y'),
                                    ('grid_longitude', 'grid_latitude')):
                for cube in self.cubes:
                    if rot_name in [i.name() for i in cube.coords(axis=ax)]:
                        cube.remove_coord(rot_name)

            try:
                _dx = thecube.attributes['um_res'].to_flt('m')
            except KeyError:
                _dx = 1.
            self.dx = Cube(_dx, units='m')
            self.dy = Cube(_dx, units='m')

        # Non-spherical coords?
        # self.horiz_cs = thecube.coord(axis='x', dim_coords=True).coord_system
        self.horiz_cs = thecube.coord_system()
        self._spherical_coords = isinstance(self.horiz_cs,
                                            (iris.coord_systems.GeogCS,
                                             iris.coord_systems.RotatedGeogCS))
        # todo: interface for spherical coordinates switch?
        # assert not self._spherical_coords,\
        #     'Only non-spherical coordinates are allowed ...'
        if self.cartesian and self._spherical_coords:
            warnings.warn('Cubes are in spherical coordinates!'
                          '\n Use `replace_lonlat_dimcoord_with_cart function`'
                          ' to change coordinates!')

    def __repr__(self):
        msg = "arke `Atmospheric Flow` containing of:\n"
        msg += "\n".join(tuple(i.name() for i in self.cubes))
        return msg

    def d_dx(self, name=None, alias=None):
        r"""
        Derivative of a cube along the x-axis
        .. math::
            \frac{\partial }{\partial x}
        """
        if name is None and isinstance(alias, str):
            v = getattr(self, alias)
        else:
            v = self.cubes.extract_strict(name)
        return cube_deriv(v, self.xcoord)

    def d_dy(self, name=None, alias=None):
        r"""
        Derivative of a cube along the y-axis
        .. math::
            \frac{\partial }{\partial y}
        """
        if name is None and isinstance(alias, str):
            v = getattr(self, alias)
        else:
            v = self.cubes.extract_strict(name)
        return cube_deriv(v, self.ycoord)

    def hgradmag(self, name=None, alias=None):
        r"""
        Magnitude of the horizontal gradient of a cube
        .. math::
            \sqrt[(\frac{\partial }{\partial y})^2
                 +(\frac{\partial }{\partial y})^2]
        """
        dvdx2 = self.d_dx(name=name, alias=alias) ** 2
        dvdy2 = self.d_dy(name=name, alias=alias) ** 2
        return (dvdx2 + dvdy2) ** 0.5

    @cached_property
    def wspd(self):
        r"""
        Calculate wind speed (magnitude)
        .. math::
            \sqrt{u^2 + v^2 + w^2}
        """
        res = 0
        for cmpnt in self.wind_cmpnt:
            res += cmpnt**2
        res = res**0.5
        res.rename('wind_speed')
        return res

    @cached_property
    def tke(self):
        r"""
        Calculate total kinetic energy
        .. math::
            0.5(u^2 + v^2 + w^2)
        """
        res = 0
        for cmpnt in self.wind_cmpnt:
            res += cmpnt**2
        res = 0.5 * res  # * self.density
        res.convert_units('m2 s-2')
        res.rename('total_kinetic_energy')
        return res

    @cached_property
    def du_dx(self):
        r"""
        Derivative of u-wind along the x-axis
        .. math::
            u_x = \frac{\partial u}{\partial x}
        """
        return cube_deriv(self.u, self.xcoord)

    @cached_property
    def du_dy(self):
        r"""
        Derivative of u-wind along the y-axis
        .. math::
            u_y = \frac{\partial u}{\partial y}
        """
        return cube_deriv(self.u, self.ycoord)

    @cached_property
    def du_dz(self):
        r"""
        Derivative of u-wind along the z-axis
        .. math::
            u_z = \frac{\partial u}{\partial z}
        """
        return cube_deriv(self.u, self.zcoord)

    @cached_property
    def dv_dx(self):
        r"""
        Derivative of v-wind along the x-axis
        .. math::
            v_x = \frac{\partial v}{\partial x}
        """
        return cube_deriv(self.v, self.xcoord)

    @cached_property
    def dv_dy(self):
        r"""
        Derivative of v-wind along the y-axis
        .. math::
            v_y = \frac{\partial v}{\partial y}
        """
        return cube_deriv(self.v, self.ycoord)

    @cached_property
    def dv_dz(self):
        r"""
        Derivative of v-wind along the z-axis
        .. math::
            v_z = \frac{\partial v}{\partial z}
        """
        return cube_deriv(self.v, self.zcoord)

    @cached_property
    def dw_dx(self):
        r"""
        Derivative of w-wind along the x-axis
        .. math::
            w_x = \frac{\partial w}{\partial x}
        """
        return cube_deriv(self.w, self.xcoord)

    @cached_property
    def dw_dy(self):
        r"""
        Derivative of w-wind along the y-axis
        .. math::
            w_y = \frac{\partial w}{\partial y}
        """
        return cube_deriv(self.w, self.ycoord)

    @cached_property
    def dw_dz(self):
        r"""
        Derivative of w-wind along the z-axis
        .. math::
            w_z = \frac{\partial w}{\partial z}
        """
        return cube_deriv(self.w, self.zcoord)

    @cached_property
    def rel_vort(self):
        r"""
        Calculate the vertical component of the vorticity vector
        .. math::
            \zeta = v_x - u_y
        """
        res = self.dv_dx - self.du_dy
        res.rename('atmosphere_relative_vorticity')
        res.convert_units('s-1')
        return res

    @cached_property
    def div_h(self):
        r"""
        Calculate the horizontal divergence
        .. math::
            D_h = u_x + v_y
        """
        res = self.du_dx + self.dv_dy
        res.rename('divergence_of_wind')
        res.convert_units('s-1')
        return res

    @cached_property
    def rel_vort_hadv(self):
        r"""
        Calculate the horizontal advection of relative vorticity
        .. math::
            \vec v\cdot \nabla \zeta
        """
        res = (self.u*cube_deriv(self.rel_vort, self.xcoord) +
               self.v*cube_deriv(self.rel_vort, self.ycoord))
        res.rename('horizontal_advection_of_atmosphere_relative_vorticity')
        res.convert_units('s-2')
        return res

    @cached_property
    def rel_vort_vadv(self):
        r"""
        Calculate the vertical advection of relative vorticity
        .. math::
            w\frac{\partial \zeta}{\partial z}
        """
        res = self.w*cube_deriv(self.rel_vort, self.zcoord)
        res.rename('vertical_advection_of_atmosphere_relative_vorticity')
        res.convert_units('s-2')
        return res

    @cached_property
    def rel_vort_stretch(self):
        r"""
        Stretching term
        .. math::
            \nabla\cdot\vec v (\zeta+f)
        """
        res = self.div_h * (self.rel_vort + self.fcor.data)
        res.rename('stretching_term_of_atmosphere_relative_vorticity_budget')
        res.convert_units('s-2')
        return res

    @cached_property
    def rel_vort_tilt(self):
        r"""
        Tilting (twisting) term
        .. math::
            \vec k \cdot \nabla w\times\frac{\partial\vec v}{\partial z} =
            \frac{\partial w}{\partial x}*\frac{\partial v}{\partial z} -
            \frac{\partial w}{\partial y}*\frac{\partial u}{\partial z}
        """
        res = self.dw_dx * self.dv_dz - self.dw_dy * self.du_dz
        res.rename('tilting_term_of_atmosphere_relative_vorticity_budget')
        res.convert_units('s-2')
        return res

    @cached_property
    def dfm_stretch(self):
        r"""
        Stretching deformation
        .. math::
            Def = u_x - v_y
        """
        res = self.du_dx - self.dv_dy
        res.rename('stretching_deformation_2d')
        res.convert_units('s-1')
        return res

    @cached_property
    def dfm_shear(self):
        r"""
        Shearing deformation
        .. math::
            Def' = u_y + v_x
        """
        res = self.du_dy + self.dv_dx
        res.rename('shearing_deformation_2d')
        res.convert_units('s-1')
        return res

    @cached_property
    def kvn(self):
        r"""
        Kinematic vorticity number

        .. math::
            W_k=\frac{||\Omega||}{||S||}=
            \frac{\sqrt{\zeta^2}}{\sqrt{D_h^2 + Def^2 + Def'^2}}
        where
        .. math::
            \zeta=v_x - u_y
            D_h = u_x + v_y
            Def = u_x - v_y
            Def' = u_y + v_x

        Reference:
            http://dx.doi.org/10.3402/tellusa.v68.29464
        """
        numerator = self.rel_vort
        denominator = (self.div_h**2 +
                       self.dfm_stretch**2 +
                       self.dfm_shear**2)**0.5
        res = numerator/denominator
        res.rename('kinematic_vorticity_number_2d')
        return res

    @cached_property
    def density(self):
        r"""
        Air density

        .. math::
            \rho = \frac{p}{R_d T}
        """
        p = self.cubes.extract_strict('air_pressure')
        rd = AuxCoord(mconst.Rd.data, units=mconst.Rd.units)
        try:
            temp = self.cubes.extract_strict('air_temperature')
        except iris.exceptions.ConstraintMismatchError:
            temp = self.temp
        res = p / (temp * rd)
        res.rename('air_density')
        res.convert_units('kg m-3')
        self.cubes.append(res)
        return res

    @cached_property
    def theta(self):
        r"""
        Air potential temperature

        If temperature is given:
        .. math::
            \theta = T (p_0/p)^{R_d/c_p}}
        """
        try:
            th = self.cubes.extract_strict('air_potential_temperature')
        except iris.exceptions.ConstraintMismatchError:
            p = self.cubes.extract_strict('air_pressure')
            temp = self.cubes.extract_strict('air_temperature')
            th = mcalc.potential_temperature(p, temp)
            th.rename('air_potential_temperature')
            th.convert_units('K')
            self.cubes.append(th)
        return th

    @cached_property
    def temp(self):
        r"""
        Air temperature

        If potential temperature is given:
        .. math::
            T = \theta (p/p_0)^{R_d/c_p}}
        """
        try:
            t = self.cubes.extract_strict('air_temperature')
        except iris.exceptions.ConstraintMismatchError:
            p0 = AuxCoord(mconst.P0.data, units=mconst.P0.units)
            kappa = mconst.kappa.data
            p = self.cubes.extract_strict('air_pressure')
            t = self.theta * (p / p0) ** kappa
            t.rename('air_temperature')
            t.convert_units('K')
            self.cubes.append(t)
        return t

    @cached_property
    def thetae(self):
        r"""
        Equivalent potential temperature

        .. math::
            \theta_e = \theta e^\frac{L_v r_s}{C_{pd} T}
        """
        try:
            th = self.cubes.extract_strict('equivalent_potential_temperature')
        except iris.exceptions.ConstraintMismatchError:
            p = self.cubes.extract_strict('air_pressure')
            temp = self.cubes.extract_strict('air_temperature')
            spechum = self.cubes.extract_strict('specific_humidity')
            mixr = spechum / ((spechum)*(-1) + 1)
            e = mcalc.vapor_pressure(p, mixr)
            dew = mcalc.dewpoint(e)
            th = mcalc.equivalent_potential_temperature(p, temp, dew)
            th.rename('equivalent_potential_temperature')
            th.convert_units('K')
            self.cubes.append(th)
        return th

    @cached_property
    def relh(self):
        r""" Relative humdity """
        try:
            rh = self.cubes.extract_strict('relative_humidity')
        except iris.exceptions.ConstraintMismatchError:
            p = self.cubes.extract_strict('air_pressure')
            temp = self.cubes.extract_strict('air_temperature')
            spechum = self.cubes.extract_strict('specific_humidity')
            rh = mcalc.specific_to_relative_humidity(p, temp, spechum)
            rh.rename('relative_humidity')
            rh.convert_units('1')
            self.cubes.append(rh)
        return rh

    @cached_property
    def specific_volume(self):
        r"""
        Air Specific Volume

        .. math::
            \alpha = \rho^{-1}
        """
        res = self.density ** (-1)
        res.rename('air_specific_volume')
        self.main_cubes.append(res)
        return res

    @cached_property
    def dp_dx(self):
        r"""
        Derivative of pressure along the x-axis
        .. math::
            p_x = \frac{\partial p}{\partial x}
        """
        return cube_deriv(self.pres, self.xcoord)

    @cached_property
    def dp_dy(self):
        r"""
        Derivative of pressure along the y-axis
        .. math::
            p_y = \frac{\partial p}{\partial y}
        """
        return cube_deriv(self.pres, self.ycoord)

    @cached_property
    def dsv_dx(self):
        r"""
        Derivative of specific volume along the x-axis
        .. math::
            \alpha_x = \frac{\partial\alpha}{\partial x}
        """
        return cube_deriv(self.specific_volume, self.xcoord)

    @cached_property
    def dsv_dy(self):
        r"""
        Derivative of specific volume along the y-axis
        .. math::
            \alpha_y = \frac{\partial\alpha}{\partial y}
        """
        return cube_deriv(self.specific_volume, self.ycoord)

    @cached_property
    def baroclinic_term(self):
        r"""
        Baroclinic term of relative vorticity budget

        .. math::
            \frac{\partial p}{\partial x}\frac{\partial\alpha}{\partial y}
            - \frac{\partial p}{\partial y}\frac{\partial\alpha}{\partial x}
        """
        res = (self.dp_dx * self.dsv_dy
               - self.dp_dy * self.dsv_dx)
        res.rename('baroclinic_term_of_atmosphere_relative_vorticity_budget')
        res.convert_units('s-1')
        return res

    @cached_property
    def brunt_vaisala_squared(self):
        r"""
        Brunt–Väisälä frequency squared

        (pressure coordinates)
        .. math::
            N^2 = -\frac{g^2 \rho}{\theta}\frac{\partial \theta}{\partial p}
        """
        dthdp = cube_deriv(self.theta, self.zcoord)
        g2 = -1 * mconst.g**2
        g2 = AuxCoord(g2.data, units=g2.units)
        res = self.density / self.theta * dthdp * g2
        res.rename('square_of_brunt_vaisala_frequency_in_air')
        res.convert_units('s-2')
        return res

    @cached_property
    def eady_growth_rate(self):
        r"""
        Eady growth rate

        (pressure coordinates)
        .. math::
            EGR = 0.31 * \frac{f}{N}\frac{\partial V}{\partial z}
                = 0.31 * \frac{-\rho g f}{N}\frac{\partial V}{\partial p}
        """
        factor = -1 * mconst.egr_factor * self.fcor * mconst.g
        dvdp = cube_deriv(self.wspd, self.zcoord)
        dvdp.data = abs(dvdp.data)
        res = (self.density * factor.data * AuxCoord(1, units=factor.units)
               * dvdp / self.brunt_vaisala_squared**0.5)
        res.rename('eady_growth_rate')
        res.convert_units('s-1')
        return res

    @cached_property
    def kinematic_frontogenesis(self):
        r"""
        2D Kinematic frontogenesis from MetPy package

        .. math::
            F=\frac{1}{2}\left|\nabla \theta\right|[D cos(2\beta)-\delta]
        """
        res = mcalc.frontogenesis(self.theta, self.u, self.v,
                                  self.dx, self.dy, dim_order='yx')
        res.rename('kinematic_frontogenesis')
        res.convert_units('K m-1 s-1')
        return res
