from math import log

from .utils import PropertiesOfWholeSystem as Prop
from fluids import friction_factor


class PropertiesOfPart:
    # properties of a part of the system
    def __init__(self, d_k=0, p_k=0, assum_increment=Prop.d * 0.1):
        # assumed fixed depth increment
        self.delta_d = assum_increment
        # assumed fixed pressure increment
        self.delta_p = Prop.p_w * 0.1
        # current depth and pressure
        self.d = d_k
        self.p = p_k
        # average depth and pressure, temperature in the k-part of the system
        self.ad = d_k + self.delta_d / 2
        self.ap = p_k + self.delta_p / 2
        self.t = self.ad * (- Prop.t_w + Prop.t_r) / Prop.d + Prop.t_w
        # Solution gas
        self.rs = 0
        # Oil formation volume factor
        self.bo = 0
        # Live oil viscosity
        self.oil_mu = 0
        # Reduced pressure and reduced temperature, dimensionless
        self.t_r, self.p_r = 0, 0
        # Gas compressibility
        self.z_factor = 0
        # Corrected volumetric flow rates,
        # where q_g for gas, q_l for liquid, q_t for total
        self.q_g, self.q_l, self.q_t = 0, 0, 0
        # Corrected mass flow rates, where w_g for gas,
        # w_l for liquid, w_t for total
        self.w_l, self.w_g, self.w_t = 0, 0, 0
        # Corrected densities, where ro_g for gas, ro_l for liquid
        self.ro_l, self.ro_g = 0, 0
        # Test variables from Appendix B. v_d - dimensionless
        # gas velocity, v_t - total fluid velocity, v_g = q_g/q_t
        self.v_t, self.v_g, self.v_d = 0, 0, 0
        # Average density, friction-loss term
        self.aver_density, self.friction_loss = 0, 0
        # Reynolds number
        self.n_re = 0
        # Moody friction factor
        self.moody_friction = 0
        # Flow type. Can be Mist, Bubble, Slug, Transition
        self.flow_type = "Bubble"
        # Boundaries for slug, mist, bubble regimes
        self.l_s, self.l_m, self.l_b = 0, 0, 0

    def dead_oil_mu(self) -> float:
        # get dead oil viscosity
        return (((- Prop.mu_d100f + Prop.mu_d210f) / 110)
                * (self.t - 100) + Prop.mu_d100f)

    def gas_oil_ratio(self) -> None:
        # Standing (1947) correlation for solution gas
        x = 0.0125 * Prop.api - 0.00091 * self.t
        self.rs = Prop.sg_g * ((self.ap / 18.2 + 1.4) * 10**x)**1.2048
        return None

    def oil_fvf(self) -> None:
        # Standing (1947) correlation for oil formation volume factor
        self.bo = (0.9759 + 0.00012 * (self.rs * (Prop.sg_g / Prop.sg_o)**0.5
                                       + 1.25 * self.t)**1.2)
        return None

    def p_bubble(self) -> float:
        # Standing’s Correlation (1947) for bubble pressure
        a = 0.00091 * self.t - 0.0125 * Prop.api
        p_bubble = 18.2 * (10**a * (self.rs / Prop.sg_g)**0.83 - 1.4)
        return p_bubble

    def gas_z_factor(self) -> None:
        # Brill and Beggs (1974) correlation for gas compressibility
        a = 1.39 * (self.t_r - 0.92)**0.5 - 0.36 * self.t_r - 0.1
        c = 0.132 - 0.32 * log(self.t_r)
        e = 9 * (self.t_r - 1)
        f = 0.3106 - 0.49 * self.t_r + 0.1824 * self.t_r**2
        d = 10**f
        b = ((0.62 - 0.23 * self.t_r) * self.p_r
             + (0.066 / (self.t_r - 0.86) - 0.037) * self.p_r**2
             + 0.32 * (self.p_r**6) / 10**e)
        self.z_factor = a + (1 - a) / e**b + c * self.p_r**d
        return None

    def live_oil_mu(self) -> None:
        # Beggs and Robinson (1975) correlation for live oil viscosity
        if self.rs > 2070:
            raise ValueError("rs should be less 2070 for "
                             "current oil viscosity correlation")
        c = 5.44 * (self.rs + 150)**-0.338
        self.oil_mu = 10.715 * (self.rs + 100)**-0.515 * self.dead_oil_mu()**c
        return None

    def tr_pr(self) -> None:
        # Reduced temperature and pressure dimensionless
        self.t_r = (self.t + 460) / Prop.t_pc
        self.p_r = self.ap / Prop.p_pc
        return None

    def volumetric_flow_rate(self) -> None:
        # corrected volumetric flow of gas, liquid, gas and liquid
        self.q_g = (3.27 * 10**-7 * self.z_factor * Prop.q_o *
                    (Prop.r - self.rs) * (self.t + 460) / self.ap)
        self.q_l = 6.49 * 10**-5 * Prop.q_o * self.bo
        self.q_t = self.q_g + self.q_l
        return None

    def velocities(self) -> None:
        # Test variables from Appendix B. v_d - dimensionless gas velocity,
        # v_t - total fluid velocity, v_g = q_g/q_t
        self.v_t = self.q_t / Prop.a
        self.v_g = self.q_g / self.q_t
        self.v_d = (self.q_g * (self.ro_l / (Prop.sigma * Prop.g))**0.25
                    / Prop.a)
        return None

    def mass_flows(self) -> None:
        # Corrected mass flow rates, where w_g for gas,
        # w_L for liquid, w_t for total
        self.w_l = (Prop.q_o * (4.05 * 10**-3 * Prop.sg_o +
                                8.85 * 10**-7 * Prop.sg_g * self.rs))
        self.w_g = 8.85 * 10**-7 * Prop.sg_g * Prop.q_o * (Prop.r - self.rs)
        self.w_t = self.w_l + self.w_g
        return None

    def corrected_densities(self) -> None:
        # Corrected densities, where ro_g for gas, ro_l for liquid
        self.ro_l = self.w_l / self.q_l
        self.ro_g = self.w_g / self.q_g
        return None

    def flow_regime(self) -> None:
        # Flow type. Can be Mist, Bubble, Slug, Transition
        self.l_b = 1.071 - (0.2218 * (self.v_t**2) / Prop.d_h)
        if self.l_b <= 0.13:
            self.l_b = 0.13
        self.l_s = 50 + 36 * (self.v_d * self.q_l) / self.q_g
        self.l_m = 75 + 84 * (self.v_g * self.q_l / self.q_g)**0.75

        if self.q_g / self.q_t < self.l_b:
            self.flow_type = "Bubble"
        elif self.q_g / self.q_t > self.l_b and self.v_d < self.l_s:
            self.flow_type = "Slug"
        elif self.l_m > self.v_t > self.l_s:
            self.flow_type = "Transition"
        elif self.v_d > self.l_m:
            self.flow_type = "Mist"
        else:
            raise "There is a problem with flow regime determination"
        return None

    def moody_friction_n_re(self, relative_roughness) -> None:
        self.n_re = 1488 * self.ro_l * Prop.d_h * self.v_t / self.oil_mu
        self.moody_friction = friction_factor(self.n_re, relative_roughness)
        return None


pass
