from math import pi


class PropertiesOfWholeSystem:
    # Constants
    g = 32.17404856  # [m/s^2] -> Acceleration of gravity, ft/sec^2
    g_c = 32.17404856  # (ft*lb)/(lb_f*s^2) -> Gravitational conversion factor
    # Fluid flow properties
    t_w = 0  # [°F] -> Temperature in the wellhead
    t_r = 0  # [°F] -> Temperature in the reservoir
    p_w = 0  # [psia] -> Pressure in the wellhead
    q_o = 0  # [STB/D] -> Oil rate
    mu_d100f = 0  # [cp] -> Dead Oil Viscosity at 100F
    mu_d210f = 0  # [cp] -> Dead Oil Viscosity at 100F
    sg_o = 0  # Oil specific gravity
    sg_g = 0  # Gas specific gravity
    a = 0  # [sq ft] -> Tubing Area
    d = 0  # [ft] -> Total depth
    r = 0  # [sqf/STB] -> Produced gas-oil ration
    d_h = 0  # [ft] -> Tubing diameter
    material = 'Carbon steel, bare'  # Material of the tube
    liq_phase = 'Oil'  # Liquid phase
    api = 0
    sigma = 0
    p_pc, t_pc = 0, 0

    @staticmethod
    def wetted_perimeter() -> float:
        # Get hydraulic pipe diameter
        p = PropertiesOfWholeSystem
        a = p.a
        return 2 * 3.14 * (a / pi) ** 0.5

    @staticmethod
    def api_grav() -> float:
        # Get API gravity
        p = PropertiesOfWholeSystem
        sg_o = p.sg_o
        if sg_o == 0:
            raise ValueError("Oil specific gravity equals zero, can't get api")
        return 141.5 / sg_o - 131.5

    @staticmethod
    def gas_pseudoprops() -> tuple[float, float]:
        # Brown et al. (1948) for gas pseudo-critical
        # properties. Only natural gas system
        p = PropertiesOfWholeSystem
        p_pc = 677 + (15 * p.sg_g) - (37.5 * p.sg_g ** 2)
        t_pc = 168 + (325 * p.sg_g) - (12.5 * p.sg_g ** 2)  # in Rankine
        return p_pc, t_pc


pass
