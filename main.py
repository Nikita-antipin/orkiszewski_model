from two_phases import TwoPhase

two_phase = TwoPhase()

# Set parameters
two_phase.t_w = 126  # [°F] -> Temperature in the wellhead
two_phase.t_r = 150  # [°F] -> Temperature in the reservoir
two_phase.p_w = 670  # [psia] -> Pressure in the wellhead
two_phase.q_o = 1850  # [STB / D] -> Oil rate
two_phase.mu_d100f = 89  # [cp] -> Dead Oil Viscosity at 100F
two_phase.mu_d210f = 8.8  # [cp] -> Dead Oil Viscosity at 210F
two_phase.sg_o = 0.942  # Oil specific gravity
two_phase.sg_g = 0.75  # Gas specific gravity
two_phase.a = 0.0488  # [sq ft] -> Tubing Area
two_phase.d = 3890  # [ft] -> Total depth
two_phase.r = 575  # [sqf / STB] -> Produced gas-oil ration
two_phase.d_h = 0.249  # [ft] -> Tubing diameter

# Get values
Depth_list, Pressure_list = two_phase.get_pressures_and_depths()
two_phase.plot()
