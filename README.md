# Orkiszewski model

This program implements Orkiszewski model and make easy to get the drop pressure and plot the graph (Depth [FT] vs Pressure [PSIA]). You are supposed to use only the single class TwoPhase() to input properties and get all important output from the model (the graph,the pressure drop).

Moody pipe relative roughness factor and friction factor are obtained from the [fluids](https://fluids.readthedocs.io/). 

## Usage

```python
from two_phases import TwoPhase

two_phase = TwoPhase()

# Set parameters
two_phase.t_w = 126  # [°F] -> Temperature in the wellhead
two_phase.t_r = 150  # [°F] -> Temperature in the reservoir
two_phase.p = 670  # [psia] -> Pressure in the wellhead (if down)
# or in the bottom (if direction is up)
two_phase.q_o = 1850  # [STB / D] -> Oil rate
two_phase.t_f = [210, 100]  # [°F] -> Temperatures, where dead oil
# viscosity were measured
two_phase.dead_mu = [8.8, 89]  # [cp] -> Dead Oil Viscosity
# at different temperature points (from t_f list)
two_phase.sg_o = 0.942  # Oil specific gravity
two_phase.sg_g = 0.75  # Gas specific gravity
two_phase.a = 0.0488  # [sq ft] -> Tubing Area
two_phase.d = 3890  # [ft] -> Total depth
two_phase.r = 575  # [sqf / STB] -> Produced gas-oil ration
two_phase.d_h = 0.249  # [ft] -> Tubing diameter
two_phase.direction = 'down'  # Direction toward the needing point
two_phase.material = 'Carbon steel, bare'  # Material of the tube
# At this version program don't support option 'Water'
two_phase.liq_phase = 'Oil'

# Get two lists of pressure and depth
depth_list, pres_list = two_phase.get_data()
# Get drop pressure and print
print(two_phase.get_pressure_drop())
# Plot graph - Depth [FT] vs Pressure [PSIA]
two_phase.plot()
```

## Plan

- [ ] Add water option


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change or add.
