from numpy import log10 as log

from fluids import friction_factor

from .utils import PropertiesOfWholeSystem as Prop
from .k_part_utils import PropertiesOfPart


class Orkiszewski:
    # This class implemented follows the paper by Orkiszewski (1967)
    # ====== Model - Extended Griffith and Wallis method by Orkiszewski =======
    @staticmethod
    def slug(p_k: PropertiesOfPart) -> tuple[float, float]:
        """Get average density and friction-loss term for slug case"""
        # The following equations are in the Appendix C, SLUG FLOW.
        # Calculation of bubble rise velocity variants
        n_re_t = 1488 * p_k.ro_l * Prop.d_h * p_k.v_t / p_k.oil_mu
        moody_friction_t = friction_factor(n_re_t, Prop.rel_rough)
        v_b3000 = ((0.546 + 8.74 * 10**-6 * n_re_t)
                   * (Prop.g * Prop.d_h) ** 0.5)
        v_b8000 = ((0.35 + 8.74 * 10**-6 * n_re_t)
                   * (Prop.g * Prop.d_h) ** 0.5)
        # Functions of bubble Reynolds number for different bubble
        # rise velocities(N_b = 1.488 * v_b * d_h * ro_l/mu_L)
        n_re_b1 = 1488 * v_b3000 * Prop.d_h * p_k.ro_l / p_k.oil_mu
        n_re_b2 = 1488 * v_b8000 * Prop.d_h * p_k.ro_l / p_k.oil_mu
        # Check what relationship we should use as
        # correct bubble rise velocity v_b
        if n_re_b1 <= 3000:
            v_b = v_b3000
        elif n_re_b2 >= 8000:
            v_b = v_b8000
        else:
            v_b_between = ((0.251 + 8.74 * 10**-6 * n_re_t)
                           * (Prop.g * Prop.d_h)**0.5)
            v_b = (0.5 * v_b_between
                   + (v_b_between**2 + 13.59 * p_k.oil_mu
                      / (p_k.ro_l * Prop.d_h**0.5))**0.5)
        # Check what function of liquid distribution coefficient (G) to use
        if p_k.v_t < 10:
            if Prop.liq_phase == 'Oil':
                g = (0.0127 * log(p_k.oil_mu + 1) / Prop.d_h**1.415 - 0.284
                     + 0.167 * log(p_k.v_t) + 0.113 * log(Prop.d_h))
            else:
                g = (0.013 * log(p_k.oil_mu) / Prop.d_h**1.38 - 0.681 + 0.232
                     * log(p_k.v_t) - 0.428 * log(Prop.d_h))
        else:
            if Prop.liq_phase == 'Oil':
                g = (0.0274 * log(p_k.oil_mu + 1) / Prop.d_h**1.371) + 0.167 \
                    + 0.569 * log(Prop.d_h) - log(p_k.v_t) \
                    * (0.01 * log(p_k.oil_mu + 1)
                       / Prop.d_h**1.571 + 0.397 + 0.63 * log(Prop.d_h))
            else:
                g = (0.045 * log(p_k.oil_mu) / Prop.d_h**0.799 - 0.709
                     - 0.162 * log(p_k.v_t) - 0.888 * log(Prop.d_h))
        # Constraint for G
        if g < - 0.065 * p_k.v_t:
            g = - 0.065 * p_k.v_t
        # Average density
        aver_density = ((p_k.w_t + p_k.ro_l * v_b * Prop.a)
                        / (p_k.q_t + v_b * Prop.a) + g * p_k.ro_l)
        # Constraint for G (C-16)
        g_cond = ((-v_b * Prop.a / (p_k.q_t + v_b * Prop.a))
                  * (1 - aver_density / p_k.ro_l))
        if p_k.v_t > 10 and g < g_cond:
            g = g_cond
            aver_density = ((p_k.w_t + p_k.ro_l * v_b * Prop.a)
                            / (p_k.q_t + v_b * Prop.a) + g * p_k.ro_l)
        # Friction-loss term
        friction_loss = (moody_friction_t * p_k.ro_l * p_k.v_t**2
                         / (2 * Prop.g_c * Prop.d_h)
                         * ((p_k.q_l + v_b * Prop.a)
                            / (p_k.q_t + v_b * Prop.a) + g))
        return aver_density, friction_loss

    @staticmethod
    def bubble(p_k: PropertiesOfPart) -> tuple[float, float]:
        """Get average density and friction-loss term for bubble case"""
        # The following equations are in the Appendix C, BUBBLE FLOW.
        # Good approximation of an average slop velocity
        # in ft/sec for bubble regime
        v_s = 0.8
        # Void fraction of gas
        fract_g = (0.5 * (1 + p_k.q_t / (v_s * Prop.a)
                          - ((1 + p_k.q_t / (v_s * Prop.a))**2
                             - 4 * p_k.q_g / (v_s * Prop.a))**0.5))
        # Average density
        aver_density = ((1 - fract_g) * p_k.ro_l
                        + fract_g * p_k.ro_g)
        v_l = p_k.q_l / (Prop.a * (1 - fract_g))
        # The Reynolds number for liquid and Moody friction
        n_re_l = 1488 * p_k.ro_l * Prop.d_h * v_l / p_k.oil_mu
        moody_friction_l = friction_factor(n_re_l, Prop.rel_rough)
        # Friction-loss term
        friction_loss_term = (moody_friction_l * v_l**2
                              / (2 * Prop.g_c * Prop.d_h))
        return aver_density, friction_loss_term

    @staticmethod
    def transition(aver_density_slug: float, fric_loss_slug: float,
                   aver_density_m: float, fric_loss_m: float,
                   p_k: PropertiesOfPart) -> tuple[float, float]:
        """Get average density and friction-loss term for transition case"""
        # The following equations are in the Appendix C, TRANSITION FLOW.
        # Slug and mist boundary limits from the Appendix B
        # Average density
        aver_density = (((p_k.l_m - p_k.v_d) / (p_k.l_m - p_k.l_s))
                        * aver_density_slug
                        + ((p_k.v_d - p_k.l_s) / (p_k.l_m - p_k.l_s))
                        * aver_density_m)
        # Friction-loss term
        frict_loss = (((p_k.l_m - p_k.v_d) / (p_k.l_m - p_k.l_s))
                      * fric_loss_slug
                      + ((p_k.v_d - p_k.l_s) / (p_k.l_m - p_k.l_s))
                      * fric_loss_m)
        return aver_density, frict_loss

    @staticmethod
    def mist(p_k: PropertiesOfPart) -> tuple[float, float]:
        """Get average density and friction-loss term for mist case"""
        # The following equations are in the appendix C, MIST FLOW.
        # new more accurate vol flow for gas for mist flow
        q_g = Prop.a * p_k.l_m * (p_k.ro_l / (Prop.g * Prop.sigma))**-0.25
        # void fraction for gas
        fract_g = 1 / (1 + p_k.q_l / q_g)
        aver_density = ((1 - fract_g) * p_k.ro_l
                        + fract_g * p_k.ro_g)
        n_w = (4.52 * 10**-7 * (p_k.v_g * p_k.oil_mu / Prop.sigma)**2
               * p_k.ro_g / p_k.ro_l)
        # Choose right relative roughness
        if n_w < 0.005:
            relative_roughness = (34 * Prop.sigma
                                  / (p_k.ro_g * Prop.d_h * p_k.v_g**2))
        else:
            relative_roughness = (174.8 * Prop.sigma * n_w**0.302
                                  / (p_k.ro_g * Prop.d_h * p_k.v_g**2))
        # The Reynolds number for gas and Moody friction
        n_re_g = 1488 * p_k.ro_g * Prop.d_h * p_k.v_g / p_k.oil_mu
        moody_friction_g = friction_factor(n_re_g, relative_roughness)
        friction_loss_term = (moody_friction_g * p_k.ro_g * p_k.v_g * p_k.v_g
                              / (2 * Prop.g_c * Prop.d_h))
        return aver_density, friction_loss_term

    @staticmethod
    def aver_density_friction_loss(p_k: PropertiesOfPart) -> None:
        """Evaluation of average density and friction-loss gradient"""
        orki = Orkiszewski
        # There is a certain way to get average density
        # and fiction-loss term for each flow type
        if p_k.flow_type == 'Slug':
            # Average density, friction-loss term
            p_k.aver_density, p_k.friction_loss = orki.slug(p_k)
        elif p_k.flow_type == 'Bubble':
            p_k.aver_density, p_k.friction_loss = orki.bubble(p_k)
        elif p_k.flow_type == 'Transition':
            aver_density_slug, fric_loss_slug = orki.slug(p_k)
            aver_density_m, fric_loss_m = orki.mist(p_k)
            p_k.aver_density, p_k.friction_loss = orki.transition(
                aver_density_slug, fric_loss_slug,
                aver_density_m, fric_loss_m, p_k)
        elif p_k.flow_type == 'Mist':
            p_k.aver_density, p_k.friction_loss = orki.mist(p_k)

    @staticmethod
    def depth_increment(p_k: PropertiesOfPart) -> float:
        """Get right depth increment"""
        return (144 * (p_k.delta_p * (1 - p_k.w_t * p_k.q_g
                                      / (4637 * Prop.a**2 * p_k.ap))
                       / (p_k.aver_density + p_k.friction_loss)))

    @staticmethod
    def pressure_drop() -> tuple[list[float], list[float]]:
        """Get two arrays of pressures and depths"""
        p_k = PropertiesOfPart(0, Prop.p) if Prop.direction == 'down'\
            else PropertiesOfPart(Prop.d, Prop.p)
        depth_list = list()
        pres_list = list()
        # Set api and d_h
        Prop.api_grav()
        Prop.d_h = Prop.wetted_perimeter() if Prop.d_h == 0 else Prop.d_h
        # Set gas pseudo-critical properties
        Prop.gas_pseudoprops()
        # Set relative roughness
        Prop.rel_roughness()
        # Surface tension, for oil, or water
        if Prop.liq_phase == "Oil":
            Prop.sigma = 0.075
        else:
            Prop.sigma = 0.159
        while 0 <= p_k.d <= Prop.d:
            # Set solution gas and pressure at bubble point
            p_k.oil_pbubble()
            p_k.gasoilratio()
            # Set pil formation volume factor
            p_k.oil_fvf()
            # Set live oil viscosity
            p_k.live_oil_mu()
            # Set reduced pressure and reduced temperature, dimensionless
            p_k.tr_pr()
            # Set gas compressibility
            p_k.gas_z_factor()
            # Set corrected volumetric flow rates,
            # where q_g for gas, q_l for liquid, q_t for total
            p_k.volumetric_flow_rate()
            # Set corrected mass flow rates,
            # where w_g for gas, w_l for liquid, w_t for total
            p_k.mass_flows()
            # Set corrected densities, where ro_g for gas, ro_l for liquid
            p_k.corrected_densities()
            # Set test variables from Appendix B. v_d - dimensionless
            # gas velocity, v_t - total fluid velocity, v_g = q_g/q_t
            p_k.velocities()
            # Set flow type. Can be Mist, Bubble, Slug, Transition
            p_k.flow_regime()
            # Set average density, friction-loss term
            Orkiszewski.aver_density_friction_loss(p_k)
            # Fixing new conditions for the second iteration.
            # New current pressure and new current depth
            new_d = p_k.d + Orkiszewski.depth_increment(p_k)
            new_p = p_k.p + p_k.delta_p
            # Check whether calculated value of depth increment
            # differ significantly (more than 5 %)
            # from assumed depth increment.
            if (abs(new_d - p_k.d) < 0.95 * abs(p_k.delta_d) or
                    abs(new_d - p_k.d) > 1.05 * abs(p_k.delta_d)):
                # Correction
                assum_increment = new_d - p_k.d
                p_k = PropertiesOfPart(p_k.d, p_k.p, assum_increment)
            else:
                depth_list.append(p_k.d)
                pres_list.append(p_k.p)
                # Next step
                p_k = PropertiesOfPart(new_d, new_p)
        return depth_list, pres_list
