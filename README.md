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
   + `d_alpha` which reports the motor rotation speed, expressed in rpm
   + `torque` which reports the motor torque, expressed in Nm
- current-speed characteristic curve, saved in the 
[`current_speed_motor_curve.csv`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/current_speed_motor_curve.csv) 
file in the [data](https://github.com/AndreaBlengino/Motor-Load-Coupling/tree/dev/data) folder.
The curve expresses the electric current absorbed by the motor as a function of its rotation speed. 
The file is structured in two columns:
   + `d_alpha` which reports the motor rotation speed, expressed in rpm
   + `current` which reports the motor torque, expressed in A
- `rotor_inertia` is the moment of inertia of the rotor. This data expressed in kgm<sup>2</sup> in the
 [`coupling_data.py`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/coupling_data.py) 
 file in the [data](https://github.com/AndreaBlengino/Motor-Load-Coupling/tree/dev/data) folder. 

The **load** is described through the following data:
- torque-position characteristic curve, saved in the 
[`load_resistance.csv`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/load_resistance.csv) 
file in the [data](https://github.com/AndreaBlengino/Motor-Load-Coupling/tree/dev/data) folder.
The curve expresses the resistant torque exerted by the load as a function of its position. 
The file is structured in two columns:
   + `beta` which reports the position of the load, expressed in °
   + `torque` which shows the resistant torque of the load, expressed in Nm
- `load_inertia` is the moment of inertia of the load. This data expressed in kgm<sup>2</sup> in the 
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

Finally it is necessary to set the initial conditions and specify some information about the time integration.
Each of the following data is saved in the 
[`coupling_data.py`](https://github.com/AndreaBlengino/Motor-Load-Coupling/blob/dev/data/coupling_data.py) file in the 
[data](https://github.com/AndreaBlengino/Motor-Load-Coupling/tree/dev/data) folder:
- `beta_0` is the initial position of the load, expressed in °
- `d_beta_0` is the initial speed of the load, expressed in rpm
- `simulation_time` is the total duration of the simulation, expressed in s
- `time_discretization` is the time discretization, expressed in s 

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

| Name                   | Variable              | Value             | Unit            |
|:-----------------------|:---------------------:|------------------:|:----------------|
| inertia of the rotor   | `rotor_inertia`       | 4x10<sup>-4</sup> | kgm<sup>2</sup> |
| inertia of the load    | `load_inertia`        | 1.2               | kgm<sup>2</sup> |
| repetitoin of the load | `load_repetition`     | `True`            |                 |
| gear ratio             | `gear_ratio`          | 80                |                 |
| gearbox efficiency     | `efficiency`          | 0.7               |                 |
| initial load position  | `beta_0`              | 0                 | °               |
| initial load speed     | `d_beta_0`            | 0                 | rpm             |
| simulation duration    | `simulation_time`     | 50                | s               |
| time discretization    | `time_discretization` | 0.001             | s               |


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
