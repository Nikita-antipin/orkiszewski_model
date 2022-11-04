from math import log

from .utils import PropertiesOfWholeSystem as Prop
from .k_part_utils import PropertiesOfPart
from fluids import roughness_Farshad, friction_factor


# здесь совсем не совпадает со статье что-то,
# там сто проц ошибка в Г и квадратные скобки это целое?


class Orkiszewski:
    # ====== Model - Extended Griffith and Wallis method by Orkiszewski =======
    #  The models implemented follows the paper by Orkiszewski (1967)
    @staticmethod
    def slug(pk: PropertiesOfPart) -> tuple[float, float]:
        # The following equations are in the Appendix C, SLUG FLOW.
        # Calculation of bubble rise velocity variants
        # to find correct one(C-7, C-8 according to the paper)
        bubble_vel3000 = ((0.546 + 8.74 * 10**-6 * pk.n_re)
                          * (Prop.g * Prop.d_h) ** 0.5)
        bubble_vel8000 = ((0.35 + 8.74 * 10**-6 * pk.n_re)
                          * (Prop.g * Prop.d_h) ** 0.5)
        # Functions of bubble Reynolds number for different bubble
        # rise velocities(N_b = 1.488 * v_b * d_h * ro_l/mu_L)
        n_re_b1 = 1488 * bubble_vel3000 * Prop.d_h * pk.ro_l / pk.oil_mu
        n_re_b2 = 1488 * bubble_vel8000 * Prop.d_h * pk.ro_l / pk.oil_mu

        # Check what relationship we should use as
        # correct bubble rise velocity v_b
        if n_re_b1 <= 3000:
            bubble_vel = bubble_vel3000
        elif n_re_b2 >= 8000:
            bubble_vel = bubble_vel8000
        else:
            bubble_vel_between = ((0.251 + 8.74 * 10**-6 * pk.n_re)
                                  * (Prop.g * Prop.d_h)**0.5)
            bubble_vel = (0.5 * bubble_vel_between
                          + (bubble_vel_between**2 + 13.59 * pk.oil_mu
                             / (pk.ro_l * Prop.d_h**0.5))**0.5)

        # Check what function of liquid distribution coefficient (G) to use
        if pk.v_t < 10:
            if Prop.liq_phase == 'Oil':
                g = (0.0127 * log(pk.oil_mu + 1) / Prop.d_h**1.415 - 0.284
                     + 0.167 * log(pk.v_t) + 0.113 * log(Prop.d_h))
            else:
                g = (0.013 * log(pk.oil_mu) / Prop.d_h**1.38 - 0.681 + 0.232
                     * log(pk.v_t) - 0.428 * log(Prop.d_h))
        else:
            if Prop.liq_phase == 'Oil':
                g = (0.0274 * log(pk.oil_mu + 1) / Prop.d_h**1.371) + 0.167 \
                    + 0.569 * log(Prop.d_h) - log(pk.v_t) \
                    * (0.01 * log(pk.oil_mu + 1)
                       / Prop.d_h**1.571 + 0.397 + 0.63 * log(Prop.d_h))
            else:
                g = (0.045 * log(pk.oil_mu) / Prop.d_h**0.799 - 0.709
                     - 0.162 * log(pk.v_t) - 0.888 * log(Prop.d_h))
        # Constraint for G
        if g < - 0.065 * pk.v_t:
            g = - 0.065 * pk.v_t

        # Average density
        aver_density = ((pk.w_t + pk.ro_l * bubble_vel * Prop.a)
                        / (pk.q_t + bubble_vel * Prop.a) + g * pk.ro_l)
        # Constraint for G
        g_cond = ((-bubble_vel * Prop.a / (pk.q_t + bubble_vel * Prop.a))
                  * (1 - aver_density / pk.ro_l))
        if pk.v_t > 10 and g < g_cond:
            g = (-bubble_vel * Prop.a / (pk.q_t + bubble_vel * Prop.a)
                 * (1 - aver_density / pk.ro_l))
            aver_density = ((pk.w_t + pk.ro_l * bubble_vel * Prop.a)
                            / (pk.q_t + bubble_vel * Prop.a) + g * pk.ro_l)
        # Friction-loss term
        friction_loss = (pk.moody_friction * pk.ro_l * pk.v_t**2
                         / (2 * Prop.g_c * Prop.d_h)
                         * ((pk.q_l + bubble_vel * Prop.a)
                            / (pk.q_t + bubble_vel * Prop.a) + g))
        return aver_density, friction_loss

    @staticmethod
    def bubble(pk: PropertiesOfPart) -> tuple[float, float]:
        # The following equations are in the Appendix C, BUBBLE FLOW.
        # Good approximation of an average slop velocity
        # in ft/sec for bubble regime
        v_s = 0.8
        # Void fraction of gas
        void_fraction_of_gas = (0.5 * (1 + pk.q_t / (v_s * Prop.a)
                                       - ((1 + pk.q_t / (v_s * Prop.a))**2
                                          - 4 * pk.q_g / (v_s * Prop.a))**0.5))
        # Average density
        aver_density = ((1 - void_fraction_of_gas) * pk.ro_l
                        + void_fraction_of_gas * pk.ro_g)
        v_l = pk.q_l / (Prop.a * (1 - void_fraction_of_gas))
        # Friction-loss term
        friction_loss_term = (pk.moody_friction * v_l * v_l
                              / (2 * Prop.g_c * Prop.d_h))
        return aver_density, friction_loss_term

    @staticmethod
    def transition(aver_density_slug: float, fric_loss_slug: float,
                   aver_density_b: float, fric_loss_b: float,
                   pk: PropertiesOfPart) -> tuple[float, float]:
        # The following equations are in the Appendix C, TRANSITION FLOW.
        # Slug and mist boundary limits from the Appendix B
        # Average density
        aver_density = (((pk.l_m - pk.v_d) / (pk.l_m - pk.l_s))
                        * aver_density_slug
                        + ((pk.v_d - pk.l_s) / (pk.l_m - pk.l_s))
                        * aver_density_b)
        # Friction-loss term
        frict_loss = (((pk.l_m - pk.v_d) / (pk.l_m - pk.l_s)) * fric_loss_slug
                      + ((pk.v_d - pk.l_s) / (pk.l_m - pk.l_s)) * fric_loss_b)
        return aver_density, frict_loss

    @staticmethod
    def mist(pk: PropertiesOfPart) -> tuple[float, float]:
        # The following equations are in the appendix C, MIST FLOW.
        # new q_g for this case
        q_g = Prop.a * pk.l_m * (pk.ro_l / (Prop.g * Prop.sigma))**-0.25
        void_fraction_of_gas = 1 / (1 + pk.q_l / q_g)
        aver_density = ((1 - void_fraction_of_gas) * pk.ro_l
                        + void_fraction_of_gas * pk.ro_g)
        n_w = (4.52 * 10**-7 * (pk.v_g * pk.oil_mu / Prop.sigma)**2
               * pk.ro_g / pk.ro_l)

        if n_w < 0.005:
            relative_roughness = (34 * Prop.sigma
                                  / (pk.ro_g * Prop.d_h * pk.v_g**2))
        else:
            relative_roughness = (174.8 * Prop.sigma * n_w**0.302
                                  / (pk.ro_g * Prop.d_h * pk.v_g**2))

        moody_friction = friction_factor(pk.n_re, relative_roughness)
        friction_loss_term = (moody_friction * pk.ro_g * pk.v_g * pk.v_g
                              / (2 * Prop.g_c * Prop.d_h))
        return aver_density, friction_loss_term

    @staticmethod
    def aver_density_friction_loss(pk: PropertiesOfPart) -> None:
        # Evaluation of average density and friction-loss gradient
        orki = Orkiszewski
        # There is a certain way to get average density
        # and fiction-loss term for each flow type
        if pk.flow_type == 'Slug':
            # Average density, friction-loss term
            aver_density, friction_loss = orki.slug(pk)
        elif pk.flow_type == 'Bubble':
            aver_density, friction_loss = orki.bubble(pk)
        elif pk.flow_type == 'Transition':
            aver_density_slug, fric_loss_slug = orki.slug(pk)
            aver_density_b, fric_loss_b = orki.bubble(pk)
            aver_density, friction_loss = orki.transition(aver_density_slug,
                                                          fric_loss_slug,
                                                          aver_density_b,
                                                          fric_loss_b, pk)
        else:
            aver_density, friction_loss = orki.mist(pk)
        pk.aver_density = aver_density
        pk.friction_loss = friction_loss
        return None

    @staticmethod
    def depth_increment(pk: PropertiesOfPart) -> float:
        return (144 * (pk.delta_p * (1 - pk.w_t * pk.q_g /
                                     (4637 * Prop.a**2 * pk.ap))
                       / (pk.aver_density + pk.friction_loss)))

    @staticmethod
    def pressure_drop() -> tuple[list[float], list[float]]:
        # Two phase pressure drops in vertical pipe
        pk = PropertiesOfPart(0, Prop.p_w)
        depth_list = list()
        pres_list = list()
        # Set api and d_h
        Prop.api = Prop.api_grav()
        Prop.d_h = Prop.wetted_perimeter() if Prop.d_h == 0 else Prop.d_h
        # Set gas pseudo-critical properties
        Prop.p_pc, Prop.t_pc = Prop.gas_pseudoprops()
        # Diameter in meters
        d = 0.3048 * 2 * (Prop.a / 3.14) ** 0.5
        # Relative roughness
        rel_rough = roughness_Farshad(Prop.material, d) / d
        # Surface tension, for oil, or water
        if Prop.liq_phase == "Oil":
            Prop.sigma = 0.075
        else:
            Prop.sigma = 0.159

        while pk.d <= Prop.d:
            # Set solution gas
            pk.gas_oil_ratio()
            # Set pil formation volume factor
            pk.oil_fvf()
            # Set live oil viscosity
            pk.live_oil_mu()
            # Set reduced pressure and reduced temperature, dimensionless
            pk.tr_pr()
            # Set gas compressibility
            pk.gas_z_factor()
            # Set corrected volumetric flow rates,
            # where q_g for gas, q_l for liquid, q_t for total
            pk.volumetric_flow_rate()
            # Set corrected mass flow rates,
            # where w_g for gas, w_l for liquid, w_t for total
            pk.mass_flows()
            # Set corrected densities, where ro_g for gas, ro_l for liquid
            pk.corrected_densities()
            # Set test variables from Appendix B. v_d - dimensionless
            # gas velocity, v_t - total fluid velocity, v_g = q_g/q_t
            pk.velocities()
            # Set flow type. Can be Mist, Bubble, Slug, Transition
            pk.flow_regime()
            # Set Reynolds number and moody friction factor
            pk.moody_friction_n_re(rel_rough)
            # Set average density, friction-loss term
            Orkiszewski.aver_density_friction_loss(pk)
            # Fixing new conditions for the second iteration.
            # New current pressure and new current depth
            new_d = pk.d + Orkiszewski.depth_increment(pk)
            new_p = pk.p + pk.delta_p
            # Check whether calculated value of depth increment
            # differ significantly from assumed depth increment.
            if (new_d - pk.d < 0.99995 * pk.delta_d or
                    new_d - pk.d > 1.00005 * pk.delta_d):
                # Correction
                assum_increment = new_d - pk.d
                pk = PropertiesOfPart(pk.d, pk.p, assum_increment)
            else:
                depth_list.append(pk.d)
                pres_list.append(pk.p)
                # Next step
                pk = PropertiesOfPart(new_d, new_p)
        return depth_list, pres_list


pass
