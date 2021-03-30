import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.interpolate import interp1d
from scipy.integrate import ode
import matplotlib.pyplot as plt
from data.coupling_data import *


class Program:

    rotor_inertia = rotor_inertia
    load_inertia = load_inertia
    load_repetition = load_repetition
    gear_ratio = gear_ratio
    efficiency = efficiency
    beta_0 = beta_0
    d_beta_0 = d_beta_0
    simulation_time = simulation_time
    step_type = step_type
    time_discretization = time_discretization
    first_step = first_step
    min_step = min_step
    max_step = max_step
    order = order

    def __init__(self):

        self.load_data()
        self.unit_conversion()
        self.interpolation_tables()
        self.array_initialization()
        self.step_type_choice()
        self.plotting()

    def load_data(self):

        # motor and load characteristics curves reading
        self.motor_curve = pd.read_csv(r'data\torque_speed_motor_curve.csv')
        self.current_curve = pd.read_csv(r'data\current_speed_motor_curve.csv')
        self.load_curve = pd.read_csv(r'data\load_resistance.csv')

    def unit_conversion(self):

        self.motor_curve['d_alpha'] = self.motor_curve['d_alpha']*2*np.pi/60       # from rpm to rad/s
        self.current_curve['d_alpha'] = self.current_curve['d_alpha']*2*np.pi/60   # from rpm to rad/s
        self.load_curve['beta'] = self.load_curve['beta']/180*np.pi                # from deg to rad
        self.beta_0 = self.beta_0/180*np.pi                                        # from deg to rad
        self.d_beta_0 = self.d_beta_0*2*np.pi/60                                   # from rpm to rad/s

    def interpolation_tables(self):

        self.motor_torque_table = interp1d(self.motor_curve['d_alpha'],
                                           self.motor_curve['torque'],
                                           fill_value = 'extrapolate')
        self.load_table = interp1d(self.load_curve['beta'],
                                   self.load_curve['torque'],
                                   fill_value = 'extrapolate')
        self.current_table = interp1d(self.current_curve['d_alpha'],
                                      self.current_curve['current'],
                                      fill_value = 'extrapolate')

    def array_initialization(self):

        self.time = [0]
        self.beta = [self.beta_0]
        self.d_beta = [self.d_beta_0]

        motor_torque = self.motor_torque_table(self.gear_ratio*self.d_beta[-1])
        resistance_torque = self.load_table(self.beta[-1])

        self.dd_beta = [(motor_torque - resistance_torque/self.gear_ratio/self.efficiency)/
                        (self.gear_ratio*self.rotor_inertia + self.load_inertia/self.gear_ratio/self.efficiency)]

        self.alpha = [self.beta[-1]*self.gear_ratio]
        self.d_alpha = [self.d_beta[-1]*self.gear_ratio]
        self.dd_alpha = [self.dd_beta[-1]*self.gear_ratio]

        self.resistance_torque = [resistance_torque]
        self.inertia_torque = [self.dd_beta[-1]*self.load_inertia]
        self.load_torque = [self.resistance_torque[-1] + self.inertia_torque[-1]]
        self.resistance_torque_to_motor = [self.load_torque[-1]/self.gear_ratio/self.efficiency]
        self.motor_torque = [motor_torque]
        self.motor_current = [self.current_table(self.d_alpha[-1])]

    def step_type_choice(self):

        if self.step_type == 'fixed':
            self.integration_fixed_step()
        elif self.step_type == 'variable':
            self.integration_variable_step()
        else:
            raise ValueError("step_type must be 'fixed' or 'variable'")

    def integration_fixed_step(self):

        for _ in tqdm(np.arange(self.time_discretization,
                                self.simulation_time + self.time_discretization,
                                self.time_discretization),
                      ncols = 100,
                      bar_format = '{desc}: {percentage:.0f}%|{bar}| {n:.0f}/{total:.0f} '
                                   '[{elapsed}<{remaining}, {rate_fmt}{postfix}]'):

            self.time.append(self.time[-1] + self.time_discretization)

            self.d_beta.append(self.d_beta[-1] + self.dd_beta[-1]*self.time_discretization)
            self.beta.append(self.beta[-1] + self.d_beta[-1]*self.time_discretization)

            self.variables_updating()

    def integration_variable_step(self):

        def ode_equation(t, y, acceleration):

            dydt = [y[1], acceleration]

            return dydt

        y_init = [self.beta[-1], self.d_beta[-1]]

        integrator = 'vode'
        params = {'method': 'bdf',
                  'nsteps': 1,
                  'first_step': self.first_step,
                  'min_step': self.min_step,
                  'max_step': self.max_step,
                  'order': self.order}

        solver = ode(ode_equation).set_integrator(integrator, **params)
        solver.set_initial_value(y_init)

        with tqdm(total = simulation_time,
                  ncols = 100,
                  bar_format = '{desc}: {percentage:.0f}%|{bar}| {n:.2f}/{total:.2f} '
                               '[{elapsed}<{remaining},{rate_fmt}{postfix}]') as progress_bar:

            while solver.t < simulation_time:

                self.time.append(solver.t)
                progress_bar.update(self.time[-1] - self.time[-2])

                solver.set_f_params(self.dd_beta[-1])
                sol = solver.integrate(simulation_time, step = True)

                self.beta.append(sol[0])
                self.d_beta.append(sol[1])

                self.variables_updating()

            progress_bar.update(simulation_time - self.time[-1])

    def variables_updating(self):

        motor_torque = self.motor_torque_table(self.gear_ratio*self.d_beta[-1])

        if self.load_repetition and self.beta[-1] > self.load_curve['beta'].max():
            resistance_torque = self.load_table(self.beta[-1]%(2*np.pi))
        else:
            resistance_torque = self.load_table(self.beta[-1])

        self.dd_beta.append((motor_torque - resistance_torque/self.gear_ratio/self.efficiency)/
                            (self.gear_ratio*self.rotor_inertia + self.load_inertia/self.gear_ratio/self.efficiency))

        self.alpha.append(self.beta[-1]*self.gear_ratio)
        self.d_alpha.append(self.d_beta[-1]*self.gear_ratio)
        self.dd_alpha.append(self.dd_beta[-1]*self.gear_ratio)

        self.resistance_torque.append(resistance_torque)
        self.inertia_torque.append(self.dd_beta[-1]*self.load_inertia)
        self.load_torque.append(self.resistance_torque[-1] + self.inertia_torque[-1])
        self.resistance_torque_to_motor.append(self.load_torque[-1]/self.gear_ratio/self.efficiency)
        self.motor_torque.append(motor_torque)
        self.motor_current.append(self.current_table(self.d_alpha[-1]))

    def plotting(self):

        axes_label_size = 14
        tick_label_size = 12
        legend_size = 14

        plt.style.use('seaborn-darkgrid')
        fig, ax = plt.subplots(3, 3, sharex = 'all')
        plt.tight_layout()
        plt.setp(ax, xlim = (self.time[0], self.time[-1]))

        ax[0, 0].plot(self.time, np.array(self.alpha)/2/np.pi)
        ax[1, 0].plot(self.time, np.array(self.d_alpha)*60/2/np.pi)
        ax[2, 0].plot(self.time, self.dd_alpha)
        ax[0, 1].plot(self.time, np.array(self.beta)/2/np.pi)
        ax[1, 1].plot(self.time, np.array(self.d_beta)*60/2/np.pi)
        ax[2, 1].plot(self.time, self.dd_beta)
        ax[0, 2].plot(self.time, self.resistance_torque, label = 'Resistant Torque')
        ax[0, 2].plot(self.time, self.inertia_torque, label = 'Inertia Torque')
        ax[0, 2].plot(self.time, self.load_torque, label = 'Load Torque')
        ax[1, 2].plot(self.time, self.resistance_torque_to_motor, label = 'Resistant Torque on Motor')
        ax[1, 2].plot(self.time, self.motor_torque, label = 'Motor Torque')
        ax[2, 2].plot(self.time, self.motor_current)

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


program = Program()
