## motor, load and gearbox data
# inertia of the rotor in kgm2
rotor_inertia = 4e-4

# inertia of the load in kgm2
load_inertia = 1.2

# repetition of the load resistant torque curve
load_repetition = True

# total gear ratio
gear_ratio = 80

# gear ratio efficiency
efficiency = 0.7


## initial conditions
# initial load position in deg
beta_0 = 0

# initial load speed in rpm
d_beta_0 = 0


## integration settings
# total simulation time in s
simulation_time = 50

# time step type: 'fixed' or 'variable'
step_type = 'fixed'

# time discretization in s, used only if step_type is 'fixed'
time_discretization = 0.001

# size of the first time step, used only if step_type is 'variable'
first_step = 1e-6

# minimum time step, used only if step_type is 'variable'
min_step = 1e-6

# maximum time step, used only if step_type is 'variable'
max_step = 1e-3

# maximum order used by the integrator, used only if step_type is 'variable'
order = 5
