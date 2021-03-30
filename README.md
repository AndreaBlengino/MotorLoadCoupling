# Script to analyze motor-load coupling dynamics

The aim of the script is to analyze the dynamics that arise during the start-up transient of a motor-load coupling.  
This coupling consists of a gearbox with a defined transmission ratio.
The calculation considers a load that exerts a resisting torque depending on its position.
Any type of electric motor with a known torque-speed characteristic can be considered as a motor. 


## Usage

### General motor-load scheme

The general scheme of the coupling is shown in the image: 

<img src="docs\scheme.png" width="830" height="300">

The position of the motor shaft is indicated by α and the load shaft one by β. 


### Input data

The **motor** is described through the following data:
- torque-speed characteristic curve, saved in the 
[`torque_speed_motor_curve.csv`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/torque_speed_motor_curve.csv) 
file in the [data](https://github.com/AndreaBlengino/Motor-Load-Coupling/tree/dev/data) folder.
The curve expresses the driving torque delivered by the motor as a function of its rotation speed. 
The file is structured in two columns:
   + `d_alpha` which reports the motor rotation speed, expressed in _rpm_
   + `torque` which reports the motor torque, expressed in _Nm_
- current-speed characteristic curve, saved in the 
[`current_speed_motor_curve.csv`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/current_speed_motor_curve.csv) 
file in the [data](https://github.com/AndreaBlengino/Motor-Load-Coupling/tree/dev/data) folder.
The curve expresses the electric current absorbed by the motor as a function of its rotation speed. 
The file is structured in two columns:
   + `d_alpha` which reports the motor rotation speed, expressed in _rpm_
   + `current` which reports the motor torque, expressed in _A_
- `rotor_inertia` is the moment of inertia of the rotor. This data expressed in _kgm<sup>2</sup>_ in the
 [`coupling_data.py`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/coupling_data.py) 
 file in the [data](https://github.com/AndreaBlengino/Motor-Load-Coupling/tree/dev/data) folder. 

The **load** is described through the following data:
- torque-position characteristic curve, saved in the 
[`load_resistance.csv`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/load_resistance.csv) 
file in the [data](https://github.com/AndreaBlengino/Motor-Load-Coupling/tree/dev/data) folder.
The curve expresses the resistant torque exerted by the load as a function of its position. 
The file is structured in two columns:
   + `beta` which reports the position of the load, expressed in _°_
   + `torque` which shows the resistant torque of the load, expressed in _Nm_
- `load_inertia` is the moment of inertia of the load. This data expressed in _kgm<sup>2</sup>_ in the 
[`coupling_data.py`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/coupling_data.py) 
file inside the [data](https://github.com/AndreaBlengino/Motor-Load-Coupling/tree/dev/data) folder.
- `load_repetition` is a boolean parameter. If, during the simulation, the value of the load position `beta` exceeds 
the range of values​expressed in the 
[`load_resistance.csv`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/load_resistance.csv) 
file you need to specify what the resisting pair should be. If `load_repetition` is `True`, the resistant torque expressed
in [`load_resistance.csv`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/load_resistance.csv) 
is repeated cyclically as the` beta` changes; if, on the other hand, it is `False`, the resistant torque does not 
repeat itself and takes place a linear extrapolation to calculate the resistant torque value. This data expressed in the 
[`coupling_data.py`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/coupling_data.py) file in the 
[data](https://github.com/AndreaBlengino/Motor-Load-Coupling/tree/dev/data) folder. 

The **gearbox** is described through the following data, each of the following data is saved in the 
[`coupling_data.py`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/coupling_data.py) file in the 
[data](https://github.com/AndreaBlengino/Motor-Load-Coupling/tree/dev/data) folder:
- `gear_ratio` is the total transmission ratio, expressed as the ratio of the rotation speed of the motor and the one of the load
- `efficiency` is the overall efficiency of the gearbox 

Then it is necessary to set the **initial conditions**.  
Each of the following data is saved in the 
[`coupling_data.py`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/coupling_data.py) file in the 
[data](https://github.com/AndreaBlengino/Motor-Load-Coupling/tree/dev/data) folder:
- `beta_0` is the initial position of the load, expressed in _°_
- `d_beta_0` is the initial speed of the load, expressed in _rpm_


Finally some parameters regarding **integration** process.  
Each of the following data is saved in the 
[`coupling_data.py`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/coupling_data.py) file in the 
[data](https://github.com/AndreaBlengino/Motor-Load-Coupling/tree/dev/data) folder:
- `simulation_time` is the total duration of the simulation, expressed in _s_
- `step_type` indicates whether to use a fixed or variable time step. Depending on its value, other parameters are set and 
two different integration methods are used:
    + if `'fixed'` the time integration is done through an [explicit Euler method](https://en.wikipedia.org/wiki/Euler_method). 
    This option is recommended in most cases.  
    In this case the user has to specify another parameter:
        * `time_discretization` is the time step, expressed in _s_
    + if `'variable'` the time integration is done through a _real-valued Variable-coefficient Ordinary Differential Equation solver 
    with a method based on backward differentiation formulas_ from 
    [`scipy.integration.ode`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html). 
    This option is recommended only if the system is governed by a [stiff equation](https://en.wikipedia.org/wiki/Stiff_equation) 
    that causes problems of numerical instability in the case of integration with a fixed time step.  
    In this case the user has to specify some other parameters:
        * `first_step` is the size of the first time step in _s_
        * `min_step` is the minimum time step in _s_
        * `max_step` is the maximum time step in _s_
        * `order`: is the maximum order used by the integrator, <= 5

#### Example

As an example, we consider a DC motor whose torque characteristic, reported in the 
[`torque_speed_motor_curve.csv`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/torque_speed_motor_curve.csv) 
file is the following:   

<img src="docs\torque_motor.png">

while the current characteristic, reported in the 
[`current_speed_motor_curve.csv`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/current_speed_motor_curve.csv) 
file is: 

<img src="docs\current_motor.png">

As an example of a load, we consider a load with the following characteristic, reported in the 
[`load_resistance.csv`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/load_resistance.csv) file: 

<img src="docs\load_torque.png">

The other data to be used as an example are shown in the table:

| Name                   | Variable              | Value             | Unit              |
|:-----------------------|:---------------------:|------------------:|:------------------|
| inertia of the rotor   | `rotor_inertia`       | 4x10<sup>-4</sup> | _kgm<sup>2</sup>_ |
| inertia of the load    | `load_inertia`        | 1.2               | _kgm<sup>2</sup>_ |
| repetitoin of the load | `load_repetition`     | `True`            |                   |
| gear ratio             | `gear_ratio`          | 80                |                   |
| gearbox efficiency     | `efficiency`          | 0.7               |                   |
| initial load position  | `beta_0`              | 0                 | _°_               |
| initial load speed     | `d_beta_0`            | 0                 | _rpm_             |
| simulation duration    | `simulation_time`     | 50                | _s_               |
| type of time step      | `step_type`           | `'fixed'`         |                   |
| time discretization    | `time_discretization` | 0.001             | _s_               |


### Output plot

The script produces a series of plots showing:
- the position, speed and acceleration of the motor
- the position, speed and acceleration of the load
- the resistant torque of the load, the inertia torque of the load and the sum of them two
- the resistant torque reported to the motor shaft, the driving torque of the motor
- the current absorbed by the motor 

#### Example

Taking the above data as an example, the calculation script will provide the following plots as a result:
 
<img src="docs\output_plot.png" width="830" height="415">


## Contributing

All contributions, bug reports, bug fixes, documentation improvements, enhancements, and ideas are welcome.


## License

[GNU General Public License v3.0](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/License.txt)
