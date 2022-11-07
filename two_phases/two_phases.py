import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from .utils import PropertiesOfWholeSystem as Prop
from .orkiszewski_model import Orkiszewski as Orkis


class TwoPhase:
    """
    This class can plot graph, output lists of data (pressure and depth lists)
    and pressure drop.
    User is supposed to use only this class.
    (Setting params is allowed only here)
    """
    # ========================== Callbacks ====================================

    @property  # t_w [°F] -> Temperature in the wellhead
    def t_w(self) -> float:
        return Prop.t_w

    @t_w.setter  # t_w [°F] -> Temperature in the wellhead
    def t_w(self, value: float):
        Prop.t_w = value

    @property  # p [psia] -> Pressure in the wellhead(if down)
    # or in the bottom(if up)
    def p(self) -> float:
        return Prop.p

    @p.setter  # p [psia] -> Pressure in the wellhead or in the bottom
    def p(self, value: float):
        Prop.p = value

    @property  # t_r [°F] -> Temperature in the reservoir
    def t_r(self) -> float:
        return Prop.t_r

    @t_r.setter  # t_r [°F] -> Temperature in the reservoir
    def t_r(self, value: float):
        Prop.t_r = value

    @property  # q_o [STB/D] -> Oil rate
    def q_o(self) -> float:
        return Prop.q_o

    @q_o.setter  # q_o [STB/D] -> Oil rate
    def q_o(self, value: float):
        Prop.q_o = value

    @property   # [cp] -> Dead Oil Viscosity
    # at different temperature points (from t_f list)
    def dead_mu(self) -> list[float]:
        return Prop.dead_mu

    @dead_mu.setter   # [cp] -> Dead Oil Viscosity
    # at different temperature points (from t_f list)
    def dead_mu(self, value: list[float]):
        Prop.dead_mu = value

    @property  # [°F] -> Temperatures, where dead oil viscosity values
    # were measured
    def t_f(self) -> list[float]:
        return Prop.t_f

    @t_f.setter  # [°F] -> Temperatures, where dead oil viscosity values
    # were measured
    def t_f(self, value: list[float]):
        Prop.t_f = value

    @property  # Oil specific gravity
    def sg_o(self) -> float:
        return Prop.sg_o

    @sg_o.setter  # Oil specific gravity
    def sg_o(self, value: float):
        Prop.sg_o = value

    @property  # Gas specific gravity
    def sg_g(self) -> float:
        return Prop.sg_g

    @sg_g.setter  # Gas specific gravity
    def sg_g(self, value: float):
        Prop.sg_g = value

    @property  # a [sq ft] -> Tubing Area
    def a(self) -> float:
        return Prop.a

    @a.setter  # a [sq ft] -> Tubing Area
    def a(self, value: float):
        Prop.a = value

    @property  # d [ft] -> Total depth
    def d(self) -> float:
        return Prop.d

    @d.setter  # d [ft] -> Total depth
    def d(self, value: float):
        Prop.d = value

    @property  # d_h [ft] -> Tubing diameter
    def d_h(self) -> float:
        return Prop.d_h

    @d_h.setter  # d_h [ft] -> Tubing diameter
    def d_h(self, value: float):
        Prop.d_h = value

    @property  # R [sqf/STB] -> Produced gas-oil ration
    def r(self) -> float:
        return Prop.r

    @r.setter  # R [sqf/STB] -> Produced gas-oil ration
    def r(self, value: float):
        Prop.r = value

    @property  # Material of the tube
    def material(self) -> str:
        return Prop.material

    @material.setter  # Material of the tube
    def material(self, value: str):
        if (value == 'Plastic coated' or value == 'Carbon steel, honed bare'
                or value == "Cr13, electropolished bare"
                or value == "Cement lining" or value == 'Carbon steel, bare'
                or value == 'Fiberglass lining' or value == 'Cr13, bare'):
            Prop.material = value
        else:
            raise ValueError("Pipe material can only be 'Plastic coated', "
                             "'Carbon steel, honed bare',"
                             "'Cr13, electropolished bare', "
                             "'Cement lining', 'Carbon steel, bare', "
                             "'Fiberglass lining', 'Cr13, bare'")

    @property  # Direction toward the needing point
    def direction(self) -> str:
        return Prop.direction

    @direction.setter  # Direction toward the needing point
    def direction(self, value: str):
        if value != 'up' and value != 'down':
            raise ValueError("Direction can only be 'up' "
                             "or 'down'")
        Prop.direction = value

    @property  # Liquid phase
    def liq_phase(self) -> str:
        return Prop.liq_phase

    @liq_phase.setter  # Liquid phase
    def liq_phase(self, value: str):
        if Prop.liq_phase != 'Oil':
            raise ValueError("Liquid phase can only be 'Water' or 'Oil. "
                             "Because for water some correlations "
                             "may not work correctly'")
        Prop.liq_phase = value
    # ================================ Output =================================

    @staticmethod
    def get_data() -> tuple[list[float], list[float]]:
        """Get lists of depth values in [FT] and pressure values in [PSIA]"""
        return Orkis.pressure_drop()

    @staticmethod
    def get_pressure_drop() -> int:
        """Get approximate pressure drop in [PSIA] at the input point"""
        depth_list, pres_list = TwoPhase.get_data()
        return pres_list[-1] - Prop.p

    @staticmethod
    def plot() -> None:
        """Plot a graph Depth - [FT] vs Pressure - [psia]"""
        depth_list, pres_list = TwoPhase.get_data()

        fig = plt.figure(figsize=(6, 9))
        ax = fig.add_subplot()

        ax.grid(linewidth=1)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(200))
        plt.gca().invert_yaxis()

        plt.xlabel("PRESSURE - PSIA")
        plt.ylabel("DEPTH - FT")
        plt.title("Calculated Pressure Drop", fontsize=15)

        ax.plot(pres_list, depth_list, color='b')
        plt.show()
