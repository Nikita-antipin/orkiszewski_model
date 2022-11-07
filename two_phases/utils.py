from numpy import pi

from fluids import roughness_Farshad


class PropertiesOfWholeSystem:
    """This class holds properties of the two-phase system as a whole"""
    # Constants
    g = 32.17404856  # [m/s^2] -> Acceleration of gravity, ft/sec^2
    g_c = 32.17404856  # (ft*lb)/(lb_f*s^2) -> Gravitational conversion factor
    # Fluid flow properties
    t_w = 0  # [°F] -> Temperature in the wellhead
    t_r = 0  # [°F] -> Temperature in the reservoir
    p = 0  # [psia] -> Pressure in the wellhead or in the bottom
    q_o = 0  # [STB/D] -> Oil rate
    t_f = [0, 0]  # [°F] -> Temperatures, where
    # dead oil viscosity were measured
    dead_mu = [0, 0]  # [cp] -> Dead Oil Viscosity at
    # different temperature points (from t_f list)
    sg_o = 0  # Oil specific gravity
    sg_g = 0  # Gas specific gravity
    a = 0  # [sq ft] -> Tubing Area
    d = 0  # [ft] -> Total depth
    r = 0  # [sqf/STB] -> Produced gas-oil ration
    d_h = 0  # [ft] -> Tubing diameter
    material = 'Carbon steel, bare'  # Material of the tube
    liq_phase = 'Oil'  # Liquid phase
    api = 0
    sigma = 0  # [lb/sec] -> Surface tension
    p_pc, t_pc = 0, 0  # [psia], [°R] -> pseudo-critical pressure, temperature
    direction = "Down"  # where is the needing point?
    rel_rough = 0  # Relative roughness

    @staticmethod
    def wetted_perimeter() -> float:
        """Get hydraulic pipe diameter"""
        p = PropertiesOfWholeSystem
        return 2 * pi * (p.a / pi) ** 0.5

    @staticmethod
    def api_grav() -> None:
        """Get API gravity"""
        p = PropertiesOfWholeSystem
        if p.sg_o == 0:
            raise ValueError("Oil specific gravity equals zero, can't get api")
        p.api = 141.5 / p.sg_o - 131.5

    @staticmethod
    def gas_pseudoprops() -> None:
        """Brown et al. (1948) for gas pseudo-critical properties"""
        # Only natural gas system
        p = PropertiesOfWholeSystem
        p.p_pc = 677 + (15 * p.sg_g) - (37.5 * p.sg_g ** 2)
        p.t_pc = 168 + (325 * p.sg_g) - (12.5 * p.sg_g ** 2)

    @staticmethod
    def rel_roughness() -> None:
        """Relative roughness in bubble and slug flow types"""
        p = PropertiesOfWholeSystem
        # Diameter in meters
        d = 0.3048 * 2 * (p.a / pi) ** 0.5
        p.rel_rough = roughness_Farshad(p.material, d) / d
