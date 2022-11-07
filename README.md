# Two-phase flow

The `two-phase` flow is a Python library that implements two-phase flows models and make easy to get the flow properties such as flow pattern, elongated bubble velocity, homogenous model properties etc. This library is structured in a way that the user can program using a simple and easy-to-use objects or in a more advanced manner can use the functions of the library directly.

The fluid properties are by default automatically obtained from the [CoolProp](https://github.com/CoolProp/CoolProp). However, you can also pass your own functions, determined experimentally or from any source you want.

The library has also some basic plot utils for some flow pattern maps.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install `two-phase` flow package.

```bash
pip install two-phase
```

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

## Roadmap

- [x] CoolProp integration
- [x] Homogeneous model
- [x] Elongated bubble models
- [x] Taitel 1980 - flow pattern map for vertical flows (being developed)
- [ ] Flow pattern map for horizontal flows
- [ ] Lockhart Martinelli model
- [ ] Alves Anular Flow model
- [ ] Taitel and Barnea model
- [ ] Drift model
- [ ] Beggs and Brill model
- [ ] Hagedorn Brown model
- [ ] Black oil model


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change or add.

When contributing, please use the [black](https://github.com/psf/black) code formatter as it formats the code to looks the  same regardless of the project you are reading.

Please make sure to update tests as appropriate.

## Credits

[@felipecastrotc](https://github.com/felipecastrotc/)

## License
[GPLv3](https://choosealicense.com/licenses/gpl-3.0/)
