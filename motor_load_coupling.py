# import
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from data.coupling_data import *

# motor and load characteristics curves reading
motor_curve = pd.read_csv(r'data\torque_speed_motor_curve.csv')
current_curve = pd.read_csv(r'data\current_speed_motor_curve.csv')
load_curve = pd.read_csv(r'data\load_resistance.csv')


# unit conversion
motor_curve['d_alpha'] = motor_curve['d_alpha']*2*np.pi/60       # from rpm to rad/s
current_curve['d_alpha'] = current_curve['d_alpha']*2*np.pi/60   # from rpm to rad/s
load_curve['beta'] = load_curve['beta']/180*np.pi                # from deg to rad

# initialize arrays
time = [0]
beta = [beta_0]
d_beta = [d_beta_0]
dd_beta = [0]
alpha = [beta[-1]*gear_ratio]
d_alpha = [d_beta[-1]*gear_ratio]
dd_alpha = [dd_beta[-1]*gear_ratio]
resistance_torque = [0]
inertia_torque = [0]
load_torque = [0]
resistance_torque_to_motor = [0]
motor_torque = [0]
motor_current = [0]

# set up interpolation tables
motor_torque_table = interp1d(motor_curve['d_alpha'], motor_curve['torque'], fill_value = 'extrapolate')
load_table = interp1d(load_curve['beta'], load_curve['torque'], fill_value = 'extrapolate')
current_table = interp1d(current_curve['d_alpha'], current_curve['current'], fill_value = 'extrapolate')

# integration
for _ in tqdm(np.arange(time_discretization, simulation_time + time_discretization, time_discretization), ncols = 100):

    time.append(time[-1] + time_discretization)

    resistance_torque_i = load_table(beta[-1])
    resistance_torque.append(resistance_torque_i)
    inertia_torque_i = dd_beta[-1]*load_inertia
    inertia_torque.append(inertia_torque_i)
    load_torque_i = resistance_torque_i + inertia_torque_i
    load_torque.append(load_torque_i)
    resistance_torque_to_motor_i = load_torque_i/gear_ratio/efficiency
    resistance_torque_to_motor.append(resistance_torque_to_motor_i)

    motor_torque_i = motor_torque_table(d_alpha[-1])
    motor_torque.append(motor_torque_i)

    dd_alpha_i = (motor_torque_i - resistance_torque_to_motor_i)/rotor_inertia
    dd_alpha.append(dd_alpha_i)
    d_alpha_i = (dd_alpha[-1] + dd_alpha[-2])*(time[-1] - time[-2])/2 + d_alpha[-1]
    d_alpha.append(d_alpha_i)
    alpha_i = (d_alpha[-1] + d_alpha[-2])*(time[-1] - time[-2])/2 + alpha[-1]
    alpha.append(alpha_i)

    motor_current_i = current_table(d_alpha_i)
    motor_current.append(motor_current_i)

    beta_i = alpha_i/gear_ratio
    beta.append(beta_i)
    d_beta_i = d_alpha_i/gear_ratio
    d_beta.append(d_beta_i)
    dd_beta_i = dd_alpha_i/gear_ratio
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
