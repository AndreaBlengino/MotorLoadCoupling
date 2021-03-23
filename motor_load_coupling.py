# import
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# motor data
I_r = 4e-4      # kgm2 - inertia of the rotor
current_curve = pd.read_csv(r'data\current_torque_motor_curve.csv')
motor_curve = pd.read_csv(r'data\torque_speed_motor_curve.csv')

# load data
I_l = 1.2       # kgm2 - inertia of the load
load_curve = pd.read_csv(r'data\load_resistance.csv')

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
motor_curve['d_alpha'] = motor_curve['d_alpha']*2*np.pi/60   # from rpm to rad/s
load_curve['beta'] = load_curve['beta']/180*np.pi            # from deg to rad

# initialize arrays
time = [0]
beta = [beta_0]
d_beta = [d_beta_0]
dd_beta = [0]
alpha = [beta[-1]*i]
d_alpha = [d_beta[-1]*i]
dd_alpha = [dd_beta[-1]*i]
resistance_torque = [0]
inertia_torque = [0]
load_torque = [0]
resistance_torque_to_motor = [0]
motor_torque = [0]
motor_current = [0]

# set up interpolation tables
motor_torque_table = interp1d(motor_curve['d_alpha'], motor_curve['torque'], fill_value = 'extrapolate')
load_table = interp1d(load_curve['beta'], load_curve['torque'], fill_value = 'extrapolate')
current_table = interp1d(current_curve['torque'], current_curve['current'], fill_value = 'extrapolate')

# integration
for _ in tqdm(np.arange(dt, T + dt, dt), ncols = 100):

    time.append(time[-1] + dt)

    resistance_torque_i = load_table(beta[-1])
    resistance_torque.append(resistance_torque_i)
    inertia_torque_i = dd_beta[-1]*I_l
    inertia_torque.append(inertia_torque_i)
    load_torque_i = resistance_torque_i + inertia_torque_i
    load_torque.append(load_torque_i)
    resistance_torque_to_motor_i = load_torque_i/i/eta
    resistance_torque_to_motor.append(resistance_torque_to_motor_i)

    motor_torque_i = motor_torque_table(d_alpha[-1])
    motor_torque.append(motor_torque_i)
    motor_current_i = current_table(motor_torque_i)
    motor_current.append(motor_current_i)

    dd_alpha_i = (motor_torque_i - resistance_torque_to_motor_i)/I_r
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
axes_label_size = 14
tick_label_size = 12
legend_size = 14

plt.style.use('seaborn-darkgrid')
fig, ax = plt.subplots(3, 3, sharex = 'all')
plt.tight_layout()
plt.setp(ax, xlim = (time[0], time[-1]))

ax[0, 0].plot(time, np.array(alpha)/2/np.pi)
ax[1, 0].plot(time, np.array(d_alpha)*60/2/np.pi)
ax[2, 0].plot(time, dd_alpha)
ax[0, 1].plot(time, np.array(beta)/2/np.pi)
ax[1, 1].plot(time, np.array(d_beta)*60/2/np.pi)
ax[2, 1].plot(time, dd_beta)
ax[0, 2].plot(time, resistance_torque, label = 'Resistance Torque')
ax[0, 2].plot(time, inertia_torque, label = 'Inertia Torque')
ax[0, 2].plot(time, load_torque, label = 'Load Torque')
ax[1, 2].plot(time, resistance_torque_to_motor, label = 'Resistance Torque on Motor')
ax[1, 2].plot(time, motor_torque, label = 'Motor Torque')
ax[2, 2].plot(time, motor_current)

ax[2, 0].set_xlabel('Time (s)')
ax[2, 1].set_xlabel('Time (s)')
ax[2, 2].set_xlabel('Time (s)')

ax[0, 0].set_ylabel(r'Motor Position, $\alpha \ (rot.)$')
ax[1, 0].set_ylabel(r'Motor Velocity, $\dot{\alpha} \ (rpm)$')
ax[2, 0].set_ylabel(r'Motor Acceleration, $\ddot{\alpha} \ (rad/s^2)$')

ax[0, 1].set_ylabel(r'Load Position, $\beta \ (rot.)$')
ax[1, 1].set_ylabel(r'Load Velocity, $\dot{\beta} \ (rpm)$')
ax[2, 1].set_ylabel(r'Load Acceleration, $\ddot{\beta} \ (rad/s^2)$')

ax[0, 2].set_ylabel(r'Torque, $T \ (Nm)$')
ax[1, 2].set_ylabel(r'Torque, $T \ (Nm)$')
ax[2, 2].set_ylabel(r'Current, $I \ (A)$')

ax[0, 2].legend(fontsize = legend_size, frameon = True)
ax[1, 2].legend(fontsize = legend_size, frameon = True)

for axi in ax.reshape(-1):
    axi.xaxis.label.set_size(axes_label_size)
    axi.yaxis.label.set_size(axes_label_size)
    for tick in axi.xaxis.get_major_ticks():
        tick.label.set_fontsize(tick_label_size)
    for tick in axi.yaxis.get_major_ticks():
        tick.label.set_fontsize(tick_label_size)

plt.show()
