import numpy as np

from numpy import exp

from scipy.optimize import fsolve  # non-linear solver

from .utils import PropertiesOfWholeSystem as Prop


class PropertiesOfPart:
    """This class holds properties of a part of the two-phase system"""
    step_variable = 0.001  # part of previous calculate pressure or depth

    def __init__(
            self, d_k=0, p_k=0, assum_increment=Prop.d * step_variable
            if Prop.direction == 'down' else -Prop.d * step_variable):
        # current depth and pressure
        self.d = d_k
        self.p = p_k
        # assumed fixed depth increment
        self.delta_d = assum_increment
        # assumed fixed pressure increment
        self.delta_p = self.p * self.step_variable if Prop.direction == 'down'\
            else -self.p * self.step_variable
        # average depth and pressure, temperature in the k-part of the system
        self.ad = d_k + self.delta_d / 2
        self.ap = p_k + self.delta_p / 2
        self.t = self.ad * (- Prop.t_w + Prop.t_r) / Prop.d + Prop.t_w
        # Solution gas
        self.r_s = 0
        # Oil formation volume factor
        self.b_o = 0
        # Live oil viscosity
        self.oil_mu = 0
        # Reduced pressure and reduced temperature, dimensionless
        self.t_r, self.p_r = 0, 0
        # Gas compressibility
        self.z = 0
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
        # Flow type. Can be Mist, Bubble, Slug, Transition
        self.flow_type = "Bubble"
        # Boundaries for slug, mist, bubble regimes
        self.l_s, self.l_m, self.l_b = 0, 0, 0
        self.p_b = 0
    # ======================= GAS - correlations ==============================

    def tr_pr(self) -> None:
        """Set reduced temperature and pressure dimensionless"""
        self.t_r = (self.t + 460) / Prop.t_pc
        self.p_r = self.ap / Prop.p_pc

    def gas_z_factor(self) -> None:
        """
        Set Gas Compressibility Factor
        * Dranchuk and Aboukassem, 1975
        """
        if 1 < self.t_r < 3 and 0.2 < self.p_r < 30:
            a1 = 0.3265
            a2 = -1.0700
            a3 = -0.5339
            a4 = 0.01569
            a5 = -0.05165
            a6 = 0.5475
            a7 = -0.7361
            a8 = 0.1844
            a9 = 0.1056
            a10 = 0.6134
            a11 = 0.7210

            def f(y):
                rho_pr, z = y
                c1 = (a1 + (a2 / self.t_r) + (a3 / (self.t_r**3))
                      + (a4 / (self.t_r**4)) + (a5 / (self.t_r**5)))
                c2 = a6 + (a7 / self.t_r) + (a8 / (self.t_r**2))
                c3 = a9 * ((a7 / self.t_r) + (a8 / (self.t_r**2)))
                c4 = (a10 * (1 + (a11 * (rho_pr**2)))
                      * ((rho_pr**2) / (self.t_r**3))
                      * (exp(-a11 * (rho_pr**2))))

                f1 = (z + (c3 * (rho_pr**5)) - (c2 * (rho_pr**2))
                      - (c1 * (rho_pr**1)) - c4 - 1)
                f2 = rho_pr - ((0.27 * self.p_r) / (z * self.t_r))
                return [f1, f2]

            roots = fsolve(f, np.array([1, 1]))  # initial guess
            self.z = roots[1]

        else:
            raise ("There is a problem with range"
                   " in correlation for gas z-factor")

    # ======================= OIL- correlations ===============================
    def oil_pbubble(self) -> None:
        """
        Set Oil Bubble-Point Pressure
        * Vazquez and Beggs, 1980
        """
        if (20 < Prop.r < 2070 and 0.56 < Prop.sg_g < 1.18
                and 16 < Prop.api < 58 and 70 < self.t < 295):
            # c1, c2, c3 coefficient from Vazquez-Beggs
            if Prop.api <= 30:
                c1 = 0.0362
                c2 = 1.0937
                c3 = 25.7240
            else:
                c1 = 0.0178
                c2 = 1.187
                c3 = 23.9310

            self.p_b = ((Prop.r / (c1 * Prop.sg_g
                                   * exp((c3 * Prop.api)
                                         / (self.t + 460))))**(1 / c2))
        else:
            raise ("There is a problem in bubble-point pressure "
                   "in the correlation")

    def dead_oil_mu(self) -> float:
        """Set dead oil viscosity"""
        return (((- Prop.dead_mu[0] + Prop.dead_mu[1])
                 / (-Prop.t_f[0] + Prop.t_f[1]))
                * (self.t - Prop.t_f[0]) + Prop.dead_mu[0])

    def gasoilratio(self) -> None:
        """
        Set gas-oil ratio
        * Standing (1947)
        """
        if self.p_b > self.ap:
            x = 0.0125 * Prop.api - 0.00091 * self.t
            self.r_s = Prop.sg_g * ((self.ap / 18.2 + 1.4) * 10**x)**1.2048
        else:
            self.r_s = Prop.r

    def oil_fvf(self) -> None:
        """
        Set Oil FVF
        * Above bubble pressure
          (Vazquez and Beggs, 1980)
        * At and less bubble pressure
          (Levitan and Murtha, 1999)
        """
        # FVF of oil at bubblepoint pressure using Levitan-Murtha
        b_o_bubble = (1 + ((0.0005 * Prop.r)
                           * ((Prop.sg_g / Prop.sg_o)**0.25))
                      + ((0.0004 * (self.t - 60)) / (Prop.sg_o * Prop.sg_g)))
        if self.ap < self.p_b:
            if Prop.api <= 30:
                # Vazquez-Beggs
                c1 = 0.0362
                c2 = 1.0937
                c3 = 25.7240
                c4 = 4.677E-4
                c5 = 1.751E-5
                c6 = -1.811E-8
            else:
                c1 = 0.0178
                c2 = 1.187
                c3 = 23.9310
                c4 = 4.670E-4
                c5 = 1.100E-5
                c6 = 1.337E-9
            r_sc = ((self.ap**c2) * c1 * Prop.sg_g
                    * exp((c3 * Prop.api) / (self.t + 460)))
            self.b_o = (1 + (c4 * r_sc) + (c5 * (self.t - 60)
                                           * (Prop.api / Prop.sg_g))
                        + (c6 * r_sc * (self.t - 60) * (Prop.api / Prop.sg_g)))
        elif self.ap == self.p_b:
            # Levitan-Murtha
            self.b_o = b_o_bubble
        else:
            # Calculate oil compressibility first using Levitan-Murtha
            coil = ((5 * Prop.r) + (17.2 * self.t) - (1180 * Prop.sg_g)
                    + (12.61 * Prop.api) - 1433) / (1E+05 * self.ap)
            # Levitan-Murtha
            self.b_o = b_o_bubble * exp(coil * (self.p_b - self.ap))

    def live_oil_mu(self) -> None:
        """
        Set live oil viscosity
        * The Beggs and Robinson (1975) Correlation
        """
        if (20 <= self.r_s <= 2070 and 0.75 <= Prop.sg_o <= 0.96
                and 0 <= self.ap <= 5250 and 70 <= self.t <= 295):
            if self.ap <= self.p_b:
                c = 5.44 * (self.r_s + 150)**-0.338
                self.oil_mu = (10.715 * (self.r_s + 100)**-0.515
                               * self.dead_oil_mu()**c)
            else:
                m = (2.6 * self.ap**1.187
                     * exp(-11.513 - 8.98 * 10**-5 * self.ap))
                self.oil_mu = self.dead_oil_mu() * (self.ap / self.p_b)**m
        else:
            raise ValueError("There is a problem with correlation "
                             "for live oil mu")
    # =========================================================================

    def volumetric_flow_rate(self) -> None:
        """Set corrected volumetric flow of gas, liquid, gas and liquid"""
        self.q_g = (3.27 * 10**-7 * self.z * Prop.q_o *
                    (Prop.r - self.r_s) * (self.t + 460) / self.ap)
        self.q_l = 6.49 * 10**-5 * Prop.q_o * self.b_o
        self.q_t = self.q_g + self.q_l

    def velocities(self) -> None:
        """Set test variables from Appendix B"""
        # v_d - dimensionless gas velocity,
        # v_t - total fluid velocity, v_g = q_g/q_t
        self.v_t = self.q_t / Prop.a
        self.v_g = self.q_g / self.q_t
        self.v_d = (self.q_g * (self.ro_l / (Prop.sigma * Prop.g))**0.25
                    / Prop.a)

    def mass_flows(self) -> None:
        """Set corrected mass flow rates"""
        # where w_g for gas,
        # w_L for liquid, w_t for total
        self.w_l = (Prop.q_o * (4.05 * 10**-3 * Prop.sg_o +
                                8.85 * 10**-7 * Prop.sg_g * self.r_s))
        self.w_g = 8.85 * 10**-7 * Prop.sg_g * Prop.q_o * (Prop.r - self.r_s)
        self.w_t = self.w_l + self.w_g

    def corrected_densities(self) -> None:
        """Set corrected densities"""
        # ro_g for gas, ro_l for liquid
        self.ro_l = self.w_l / self.q_l
        self.ro_g = self.w_g / self.q_g

    def flow_regime(self) -> None:
        """Flow type. Can be Mist, Bubble, Slug, Transition"""
        self.l_b = 1.071 - 0.2218 * self.v_t**2 / Prop.d_h
        if self.l_b <= 0.13:
            self.l_b = 0.13
        self.l_s = 50 + 36 * self.v_d * self.q_l / self.q_g
        self.l_m = 75 + 84 * (self.v_g * self.q_l / self.q_g)**0.75

        if self.q_g / self.q_t < self.l_b:
            self.flow_type = "Bubble"
        elif self.q_g / self.q_t > self.l_b and self.v_d < self.l_s:
            self.flow_type = "Slug"
        elif self.l_m > self.v_d > self.l_s:
            self.flow_type = "Transition"
        elif self.v_d > self.l_m:
            self.flow_type = "Mist"
        else:
            raise "There is a problem with flow regime determination"
