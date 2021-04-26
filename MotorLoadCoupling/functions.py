import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import ode
from tqdm import tqdm
import matplotlib.pyplot as plt


def load_data(parameters):

    parameters['motor_torque_curve'] = pd.read_csv(parameters['motor_torque_curve_file'])
    parameters['motor_current_curve'] = pd.read_csv(parameters['motor_current_curve_file'])
    parameters['load_torque_curve'] = pd.read_csv(parameters['load_torque_curve_file'])

def convert_units(parameters):

    parameters['motor_torque_curve']['d_alpha'] = parameters['motor_torque_curve']['d_alpha']*2*np.pi/60
    parameters['motor_current_curve']['d_alpha'] = parameters['motor_current_curve']['d_alpha']*2*np.pi/60
    parameters['load_torque_curve']['beta'] = parameters['load_torque_curve']['beta']/180*np.pi
    parameters['beta_0'] = parameters['beta_0']/180*np.pi
    parameters['d_beta_0'] = parameters['d_beta_0']*2*np.pi/60

def calculate_interpolation_tables(parameters):

    parameters['motor_torque_table'] = interp1d(parameters['motor_torque_curve']['d_alpha'],
                                                parameters['motor_torque_curve']['torque'],
                                                fill_value = 'extrapolate')
    parameters['motor_current_table'] = interp1d(parameters['motor_current_curve']['d_alpha'],
                                                 parameters['motor_current_curve']['current'],
                                                 fill_value = 'extrapolate')
    parameters['load_torque_table'] = interp1d(parameters['load_torque_curve']['beta'],
                                               parameters['load_torque_curve']['torque'],
                                               fill_value = 'extrapolate')

def initialize_array(parameters, variables):

    variables['time'] = [0]
    variables['beta'] = [parameters['beta_0']]
    variables['d_beta'] = [parameters['d_beta_0']]

    variables['motor_torque'] = [parameters['motor_torque_table'](parameters['gear_ratio']*variables['d_beta'][-1]).item()]
    variables['load_torque'] = [parameters['load_torque_table'](variables['beta'][-1]).item()]

    variables['dd_beta'] = [(variables['motor_torque'][-1] - variables['load_torque'][-1]/
                             parameters['gear_ratio']/parameters['efficiency'])/
                            (parameters['gear_ratio']*parameters['rotor_inertia'] + parameters['load_inertia']/
                             parameters['gear_ratio']/parameters['efficiency'])]

    variables['alpha'] = [variables['beta'][-1]*parameters['gear_ratio']]
    variables['d_alpha'] = [variables['d_beta'][-1]*parameters['gear_ratio']]
    variables['dd_alpha'] = [variables['dd_beta'][-1]*parameters['gear_ratio']]

    variables['inertia_torque'] = [variables['dd_beta'][-1]*parameters['load_inertia']]
    variables['resistant_torque'] = [variables['load_torque'][-1] + variables['inertia_torque'][-1]]
    variables['resistant_torque_to_motor'] = [variables['resistant_torque'][-1]/
                                              parameters['gear_ratio']/parameters['efficiency']]

    variables['motor_current'] = [parameters['motor_current_table'](variables['d_alpha'][-1]).item()]

def fixed_step_time_integration(parameters, variables):

    for _ in tqdm(np.arange(parameters['time_discretization'],
                            parameters['simulation_time']+ parameters['time_discretization'],
                            parameters['time_discretization']),
                  ncols = 100,
                  bar_format = '{desc}: {percentage:.0f}%|{bar}| {n:.0f}/{total:.0f} '
                               '[{elapsed}<{remaining}, {rate_fmt}{postfix}]'):

        variables['time'].append(variables['time'][-1] + parameters['time_discretization'])

        variables['d_beta'].append(variables['d_beta'][-1] + variables['dd_beta'][-1]*parameters['time_discretization'])
        variables['beta'].append(variables['beta'][-1] + variables['d_beta'][-1]*parameters['time_discretization'])

        variables_updating(parameters, variables)

def variable_step_time_integration(parameters, variables):

    def ode_equation(t, y, acceleration):
        """
        Defines ODE equation. Returns dydt which is a list of the derivatives of y:
        the first element is the speed, the second one the acceleration.
        """

        dydt = [y[1], acceleration]

        return dydt

    y_init = [variables['beta'][-1], variables['d_beta'][-1]]

    integrator = 'vode'
    params = {'method': 'bdf',
              'nsteps': 1,
              'first_step': parameters['first_step'],
              'min_step': parameters['min_step'],
              'max_step': parameters['max_step'],
              'order': parameters['order']}

    solver = ode(ode_equation).set_integrator(integrator, **params)
    solver.set_initial_value(y_init)

    with tqdm(total = parameters['simulation_time'],
              ncols = 100,
              bar_format = '{desc}: {percentage:.0f}%|{bar}| {n:.2f}/{total:.2f} '
                           '[{elapsed}<{remaining},{rate_fmt}{postfix}]') as progress_bar:

        while solver.t < parameters['simulation_time']:

            variables['time'].append(solver.t)
            progress_bar.update(variables['time'][-1] - variables['time'][-2])

            solver.set_f_params(variables['dd_beta'][-1])
            sol = solver.integrate(parameters['simulation_time'], step = True)

            variables['beta'].append(sol[0])
            variables['d_beta'].append(sol[1])

            variables_updating(parameters, variables)

        progress_bar.update(parameters['simulation_time'] - variables['time'][-1])

def variables_updating(parameters, variables):

    variables['motor_torque'].append(parameters['motor_torque_table'](parameters['gear_ratio']*variables['d_beta'][-1]))

    if parameters['load_repetition'] and variables['beta'][-1] > parameters['load_torque_curve']['beta'].max():
        variables['load_torque'].append(parameters['load_torque_table'](variables['beta'][-1]%(2*np.pi)))
    else:
        variables['load_torque'].append(parameters['load_torque_table'](variables['beta'][-1]))

    variables['dd_beta'].append((variables['motor_torque'][-1] - variables['load_torque'][-1]/
                                 parameters['gear_ratio']/parameters['efficiency'])/
                                (parameters['gear_ratio']*parameters['rotor_inertia'] +
                                 parameters['load_inertia']/parameters['gear_ratio']/parameters['efficiency']))

    variables['alpha'].append(variables['beta'][-1]*parameters['gear_ratio'])
    variables['d_alpha'].append(variables['d_beta'][-1]*parameters['gear_ratio'])
    variables['dd_alpha'].append(variables['dd_beta'][-1]*parameters['gear_ratio'])

    variables['inertia_torque'].append(variables['dd_beta'][-1]*parameters['load_inertia'])
    variables['resistant_torque'].append(variables['load_torque'][-1] + variables['inertia_torque'][-1])
    variables['resistant_torque_to_motor'].append(variables['resistant_torque'][-1]/
                                                  parameters['gear_ratio']/parameters['efficiency'])
    variables['motor_current'].append(parameters['motor_current_table'](variables['d_alpha'][-1]))

def plot_variables(variables):

    axes_label_size = 14
    tick_label_size = 12
    legend_size = 14

    plt.style.use('seaborn-darkgrid')
    fig, ax = plt.subplots(3, 3, sharex = 'all')
    plt.tight_layout()
    plt.setp(ax, xlim = (variables['time'][0], variables['time'][-1]))

    ax[0, 0].plot(variables['time'], np.array(variables['alpha'])/2/np.pi)
    ax[1, 0].plot(variables['time'], np.array(variables['d_alpha'])*60/2/np.pi)
    ax[2, 0].plot(variables['time'], variables['dd_alpha'])
    ax[0, 1].plot(variables['time'], np.array(variables['beta'])/2/np.pi)
    ax[1, 1].plot(variables['time'], np.array(variables['d_beta'])*60/2/np.pi)
    ax[2, 1].plot(variables['time'], variables['dd_beta'])
    ax[0, 2].plot(variables['time'], variables['load_torque'], label = 'Load Torque')
    ax[0, 2].plot(variables['time'], variables['inertia_torque'], label = 'Inertia Torque')
    ax[0, 2].plot(variables['time'], variables['resistant_torque'], label = 'Resistant Torque')
    ax[1, 2].plot(variables['time'], variables['resistant_torque_to_motor'], label = 'Resistant Torque on Motor')
    ax[1, 2].plot(variables['time'], variables['motor_torque'], label = 'Motor Torque')
    ax[2, 2].plot(variables['time'], variables['motor_current'])

    ax[2, 0].set_xlabel('Time (s)')
    ax[2, 1].set_xlabel('Time (s)')
    ax[2, 2].set_xlabel('Time (s)')

    ax[0, 0].set_ylabel(r'Motor Position, $\alpha \ (rot.)$')
    ax[1, 0].set_ylabel(r'Motor Speed, $\dot{\alpha} \ (rpm)$')
    ax[2, 0].set_ylabel(r'Motor Acceleration, $\ddot{\alpha} \ (rad/s^2)$')

    ax[0, 1].set_ylabel(r'Load Position, $\beta \ (rot.)$')
    ax[1, 1].set_ylabel(r'Load Speed, $\dot{\beta} \ (rpm)$')
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
