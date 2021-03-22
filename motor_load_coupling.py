# import
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# motor data
n_0 = 3000      # rpm  - no load veloctiy
T_s = 30        # Nmm  - stall torque
I_r = 4e-4      # kgm2 - inertia of the rotor
I_0 = 0.043     # A    - no load current
I_s = 0.87      # A    - stall current

# load data
I_l = 1.2       # kgm2 - inertia of the load
df = pd.read_csv(r'data\load_resistance.csv')

# gearbox data
i = 80          #      - total gear ratio
eta = 0.7       #      - gearbox efficiency

# initial conditions
beta_0 = 0     # rad    - initial load position
d_beta_0 = 0   # rad/s  - initial load velocity

# simulation time
T = 50           # s   - simulation time
dt = 0.001       # s   - time discretization


# unit conversion
n_0 = n_0*2*np.pi/60                # from rpm to rad/s
T_s = T_s/1000                      # from Nmm to Nm
df['beta'] = df['beta']/180*np.pi   # from deg to rad

# initialize arrays
time = [0]
beta = [beta_0]
d_beta = [d_beta_0]
dd_beta = [0]
alpha = [beta[-1]*i]
d_alpha = [d_beta[-1]*i]
dd_alpha = [dd_beta[-1]*i]

# set up interpolation tables
load_table = interp1d(df['beta'], df['torque'], fill_value = 'extrapolate')
motor_torque_table = interp1d([0, n_0], [T_s, 0], fill_value = 'extrapolate')


# integration
for _ in tqdm(np.arange(dt, T + dt, dt), ncols = 100):

    time.append(time[-1] + dt)

    resistance_torque = load_table(beta[-1])
    inertia_torque = dd_beta[-1]*I_l
    load_torque = resistance_torque + inertia_torque

    resistance_torque_to_motor = load_torque/i/eta
    motor_torque = motor_torque_table(d_alpha[-1])
    dd_alpha_i = (motor_torque - resistance_torque_to_motor)/I_r
    dd_alpha.append(dd_alpha_i)

    d_alpha_i = (dd_alpha[-1] + dd_alpha[-2])*(time[-1] - time[-2])/2 + d_alpha[-1]
    d_alpha.append(d_alpha_i)

    alpha_i = (d_alpha[-1] + d_alpha[-2])*(time[-1] - time[-2])/2 + alpha[-1]
    alpha.append(alpha_i)

    beta_i = alpha_i/i
    beta.append(beta_i)

    d_beta_i = d_alpha_i/i
    d_beta.append(d_beta_i)

    dd_beta_i = dd_alpha_i/i
    dd_beta.append(dd_beta_i)


# plotting
fig, ax = plt.subplots(3, 2)

ax[0, 0].plot(time, alpha)
ax[1, 0].plot(time, d_alpha)
ax[2, 0].plot(time, dd_alpha)
ax[0, 1].plot(time, beta)
ax[1, 1].plot(time, d_beta)
ax[2, 1].plot(time, dd_beta)

plt.show()
