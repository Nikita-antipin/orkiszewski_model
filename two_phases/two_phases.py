import matplotlib.pyplot as plt
from .utils import PropertiesOfWholeSystem as Prop
from .orkiszewski_model import Orkiszewski as Orkis

# здесь совсем не совпадает со статье что-то,
# там сто проц ошибка в Г и квадратные скобки это целое?


class TwoPhase:
    # ======================== Callbacks ================================

    @property  # t_w [°F] -> Temperature in the wellhead
    def t_w(self) -> float:
        return Prop.t_w

    @t_w.setter  # t_w [°F] -> Temperature in the wellhead
    def t_w(self, value: float):
        Prop.t_w = value

    @property  # p_w [psia] -> Pressure in the wellhead
    def p_w(self) -> float:
        return Prop.p_w

    @p_w.setter  # p_w [psia] -> Pressure in the wellhead
    def p_w(self, value: float):
        Prop.p_w = value

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

    @property   # mu_d100F [cp] -> Dead Oil Viscosity at 100F
    def mu_d100f(self) -> float:
        return Prop.mu_d100f

    @mu_d100f.setter   # mu_d100F [cp] -> Dead Oil Viscosity at 100F
    def mu_d100f(self, value: float):
        Prop.mu_d100f = value

    @property  # mu_d210F [cp] -> Dead Oil Viscosity at 210F
    def mu_d210f(self) -> float:
        return Prop.mu_d210f

    @mu_d210f.setter  # mu_d210F [cp] -> Dead Oil Viscosity at 210F
    def mu_d210f(self, value: float):
        Prop.mu_d210f = value

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
        if value != ('Plastic coated' or value != 'Carbon steel, honed bare' or
                     value != "Cr13, electropolished bare" or
                     value != "Cement lining" or
                     value != 'Carbon steel, bare' or 'Fiberglass lining' or
                     'Cr13, bare'):
            raise ValueError("Pipe material can only be 'Plastic coated', "
                             "'Carbon steel honed bare',"
                             "'Cr13, electropolished bare', "
                             "'Cement lining', 'Carbon steel, bare', "
                             "'Fiberglass lining', 'Cr13, bare'")
        Prop.material = value

    @property  # Liquid phase
    def liq_phase(self) -> str:
        return Prop.liq_phase

    @liq_phase.setter  # Liquid phase
    def liq_phase(self, value: str):
        if Prop.liq_phase != 'Water' or Prop.liq_phase != 'Oil':
            raise ValueError("Liquid phase can only be 'Water' or 'Oil'")
        Prop.liq_phase = value
    # ================================ Output =================================

    @staticmethod
    def get_pressures_and_depths() -> tuple[list[float], list[float]]:
        """
        Get list of pressure values and depth values
        """
        return Orkis.pressure_drop()

    @staticmethod
    def plot() -> None:
        """
        Plot a graph Depth - [FT] vs Pressure - [psia]
        """
        tp = TwoPhase
        plt.gca().invert_yaxis()
        depth_list, pres_list = tp.get_pressures_and_depths()
        plt.scatter(pres_list, depth_list)
        plt.show()
        return None


pass
