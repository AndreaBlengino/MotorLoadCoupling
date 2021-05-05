# Motor-Gearbod-Load Coupling Dynamics Analysis

**MotorLoadCoupling** is a python library for analyzing the start-up
transient dynamics of a motor-gearbox-load coupling.  
The calculation considers a load that exerts a resistant torque
depending on its position. Any type of electric motor with a known
torque-speed characteristic can be considered as a motor.


## Index

- [Usage](#usage)
  + [General motor-load scheme](#general-motor-load-scheme)
  + [Quick Start](#quick-start)
  + [Input data](#input-data)
  + [Results plot](#results-plot)
  + [Example](#example)
- [Requirements](#requirements)
- [Python versions](#python-versions)
- [Contributing](#contributing)
- [License](#license)
  
## Usage

### General motor-load scheme

The general scheme of the coupling is shown in the image: 

<img src="docs\conceptual_scheme.png" width="830" height="270">

The position of the motor shaft is indicated by α and the load shaft one
by β.  
The integration scheme of the coupling is shown in the image:

 <img src="docs\mathematical_scheme.png" width="830" height="340">


### Quick Start


```python
from MotorLoadCoupling import Coupling
```

Initialize class

```python
coupling = Coupling()
```


### Input data

`set_characteristics_curves_files` specifies the paths to the files
where motor and load characteristics curves are stored. These files must
be in .csv format.

```python
coupling.set_characteristics_curves_files(motor_torque_curve_file = r'data/torque_speed_motor_curve.csv',
                                          motor_current_curve_file = r'data/current_speed_motor_curve.csv',
                                          load_torque_curve_file = r'data/torque_position_load_curve.csv')
```

- `motor_torque_curve_file` points to a file where the **motor's 
  torque-speed characteristic curve** is stored. The curve expresses the
  driving torque delivered by the motor as a function of its rotation
  speed. The file is structured in two columns:
    + `d_alpha` which reports the motor rotation speed, expressed in 
    _rpm_
    + `torque` which reports the motor torque, expressed in _Nm_

- `motor_current_curve_file` points to a file where the **motor's
  current-speed characteristic curve** is stored. The curve
  expresses the electric current absorbed by the motor as a function of
  its rotation speed. The file is structured in two columns:
   + `d_alpha` which reports the motor rotation speed, expressed in 
   _rpm_
   + `current` which reports the motor torque, expressed in _A_
    
- `load_torque_curve_file` points to a file where the **load's
  torque-position characteristic curve** is stored. The curve expresses
  the resistant torque exerted by the load as a function of its
  position. The file is structured in two columns:
   + `beta` which reports the position of the load, expressed in _°_
   + `torque` which shows the resistant torque of the load, expressed in
    _Nm_

`set_inertia` specifies motor's rotor and load's **moments of inertia**.
Values are in _kgm<sup>2</sup>_.

```python
coupling.set_inertia(rotor_inertia = 4e-4, load_inertia = 1.2)
```

`set_gearbox` specifies **gearbox** parameters as efficiency and the
global gear ratio (expressed as the ratio of the rotation speed of the
motor and the one of the load).  
`load_repetition` is a boolean parameter. If, during the simulation, the
value of the load position `beta` exceeds the range of values​expressed
in the `load_torque_curve_file` characteristic curve, then it is
necessary to specify what the resistant torque value should be. If
`load_repetition` is `True`, then the resistant torque curve expressed
in `load_torque_curve_file` is repeated cyclically as the `beta`
changes; if, on the other hand, it is `False`, then the resistant torque
does not repeat itself and the resistant torque is calculated by linear
extrapolation.

```python 
coupling.set_gearbox(load_repetition = True, gear_ratio = 80, efficiency = 0.7)
```

The load's **initial position and velocity** are specified with
`set_initial_conditions`. `beta_0` is expressed in _°_ and `d_beta_0` in
_rpm_.

```python 
coupling.set_initial_conditions(beta_0 = 0, d_beta_0 = 0)
```

Then it is necessary to specify the **total duration** of the simulation
in _s_ and the **step type**.  
`step_type` can be only `'fixed'` or `'variable'`. If it is`'fixed'`,
the time integration is done through an
[explicit Euler method](https://en.wikipedia.org/wiki/Euler_method);
this option is recommended in most cases. If it is `'variable'`, the
time integration is done through a _real-valued Variable-coefficient
Ordinary Differential Equation solver with a method based on backward
differentiation formulas_ from
[`scipy.integration.ode`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html#scipy-integrate-ode).
This option is recommended only if the system is governed by a
[stiff equation](https://en.wikipedia.org/wiki/Stiff_equation) that
causes problems of numerical instability in the case of integration with
a fixed time step.

```python 
coupling.set_time(simulation_time = 50, step_type = 'fixed')
```

`set_integrator` sets the integrator parameters. The required parameters
depend on `step_type` value:

- if `step_type` is `'fixed'`, then it is necessary to specify the
  discretization time step, in _s_.
  
  ```python
  coupling.set_integrator(time_discretization = 0.01)
  ```
  
- if `step_type` is `'variable'`, then it is necessary to specify the
  first, the minimum and the maximum time step and the maximum order of
  the integrator (<=5).
  
  ```python 
  coupling.set_integrator(first_step = 1e-6, min_step = 1e-6, max_step = 1e-3, order = 5)
  ```
  
To start simulation:

```python
coupling.start()
```

All simulation input parameters are stored in the `coupling.parameters`
dictionary:

```python
for param, value in coupling.parameters.items():
    print(f'{param} = {value}')
```

```python
motor_torque_curve_file = data/torque_speed_motor_curve.csv
motor_current_curve_file = data/current_speed_motor_curve.csv
load_torque_curve_file = data/torque_position_load_curve.csv
rotor_inertia = 0.0004
load_inertia = 1.2
load_repetition = True
gear_ratio = 80
...
```

All calculated time-dependent variables are stored in the
`coupling.variables` dictionary:

```python
print(coupling.variables['beta'])
```

```python
[0.0, 5.614973262032086e-05, 0.00016836891265989472, ..., 153.163131795369, 153.19863939816887, 153.23413900522692]
```

Optionally, it is possible to plot the time-dependent variables with the
method:

```python
coupling.plotting()
```


### Results plot

The script produces a series of plots showing:
- the position, speed and acceleration of the motor
- the position, speed and acceleration of the load
- the load torque, the inertia torque of the load and the sum of them 
  two
- the resistant torque reported to the motor shaft, the driving torque 
  of the motor
- the current absorbed by the motor 


### Example

Consider the above input data and the following motor and load
characteristics curves:

- DC **motor's torque-speed characteristic curve** reported in
  [`data\torque_speed_motor_curve.csv`](https://github.com/AndreaBlengino/MotorLoadCoupling/blob/master/data/torque_speed_motor_curve.csv):

<img src="docs\torque_speed_motor_curve.png">

- DC **motor's current-speed characteristic curve** reported in
  [`data\current_speed_motor_curve.csv`](https://github.com/AndreaBlengino/MotorLoadCoupling/blob/master/data/current_speed_motor_curve.csv):

<img src="docs\current_speed_motor_curve.png">

- **load's torque-position characteristic curve** reported in
  [`data\torque_position_load_curve.csv`](https://github.com/AndreaBlengino/MotorLoadCoupling/blob/master/data/torque_position_load_curve.csv):

<img src="docs\torque_position_load_curve.png">

The simulation parameters used for an example and already described
above are:

| Name                   | Variable              | Value             | Unit              |
|:-----------------------|:---------------------:|------------------:|:------------------|
| inertia of the rotor   | `rotor_inertia`       | 4×10<sup>-4</sup> | _kgm<sup>2</sup>_ |
| inertia of the load    | `load_inertia`        | 1.2               | _kgm<sup>2</sup>_ |
| repetitoin of the load | `load_repetition`     | `True`            |                   |
| gear ratio             | `gear_ratio`          | 80                |                   |
| gearbox efficiency     | `efficiency`          | 0.7               |                   |
| initial load position  | `beta_0`              | 0                 | _°_               |
| initial load speed     | `d_beta_0`            | 0                 | _rpm_             |
| simulation duration    | `simulation_time`     | 50                | _s_               |
| type of time step      | `step_type`           | `'fixed'`         |                   |
| time discretization    | `time_discretization` | 0.01              | _s_               |


Taking the above data as an example, the time integration provides the 
following plots as a result:
 
<img src="docs\results_plot.png" width="830" height="415">


## Requirements

- [matplotlib](https://matplotlib.org) >= 3.0.1
- [numpy](https://numpy.org) >= 1.19.5
- [pandas](https://pandas.pydata.org) >= 0.25.3
- [scipy](https://scipy.org) >= 1.5.4
- [tqdm](https://tqdm.github.io) >= 4.59.0


## Python versions

MotorLoadCoupling runs normal for python versions 3.6+.


## Contributing

All contributions, bug reports, bug fixes, documentation improvements, 
enhancements, and ideas are welcome.


## License

[GNU General Public License v3.0](https://github.com/AndreaBlengino/MotorLoadCoupling/blob/master/LICENSE)
