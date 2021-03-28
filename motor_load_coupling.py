import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.interpolate import interp1d
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
    time_discretization = time_discretization

    def __init__(self):

        self.load_data()
        self.unit_conversion()
        self.interpolation_tables()
        self.dd_beta_iterative_correction()
        self.array_initialization()
        self.integration()
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
        self.alpha = [self.beta[-1]*self.gear_ratio]
        self.d_alpha = [self.d_beta[-1]*self.gear_ratio]
        self.dd_alpha = [self.dd_beta[-1]*self.gear_ratio]
        self.resistance_torque = [self.load_table(self.beta[-1])]
        self.inertia_torque = [self.dd_beta[-1]*self.load_inertia]
        self.load_torque = [self.resistance_torque[-1] + self.inertia_torque[-1]]
        self.resistance_torque_to_motor = [self.load_torque[-1]/self.gear_ratio/self.efficiency]
        self.motor_torque = [self.motor_torque_table(self.d_alpha[-1])]
        self.motor_current = [self.current_table(self.d_alpha[-1])]

    def dd_beta_iterative_correction(self):

        # first attempt value
        self.dd_beta = [0]

        for _ in range(20):

            self.array_initialization()

            self.dd_alpha.append((self.motor_torque[-1] - self.resistance_torque_to_motor[-1])/self.rotor_inertia)
            self.d_alpha.append((self.dd_alpha[-1] + self.dd_alpha[-2])*self.time_discretization/2 + self.d_alpha[-1])
            self.d_beta.append(self.d_alpha[-1]/self.gear_ratio)

            self.dd_beta_0 = (self.d_beta[-1] - self.d_beta[-2])/self.time_discretization
            self.dd_beta = [self.dd_beta_0]

    def integration(self):

        for _ in tqdm(np.arange(self.time_discretization,
                                self.simulation_time + self.time_discretization,
                                self.time_discretization), ncols = 100):

            self.time.append(self.time[-1] + self.time_discretization)

            # calculation of the resistant torque to the motor by adding the load resistance
            # by its characteristic curve and the load inertia, all divided by the total gear ratio
            if self.load_repetition and self.beta[-1] > self.load_curve['beta'].max():
                resistance_torque_i = self.load_table(self.beta[-1] % (2*np.pi))
            else:
                resistance_torque_i = self.load_table(self.beta[-1])
            self.resistance_torque.append(resistance_torque_i)
            inertia_torque_i = self.dd_beta[-1]*self.load_inertia
            self.inertia_torque.append(inertia_torque_i)
            load_torque_i = resistance_torque_i + inertia_torque_i
            self.load_torque.append(load_torque_i)
            resistance_torque_to_motor_i = load_torque_i/self.gear_ratio/self.efficiency
            self.resistance_torque_to_motor.append(resistance_torque_to_motor_i)

            # calculation of the motor torque by its characteristic curve
            motor_torque_i = self.motor_torque_table(self.d_alpha[-1])
            self.motor_torque.append(motor_torque_i)

            # calculation of the motor acceleration by torque equilibrium equation
            # and motor speed and position by integration with trapezoidal rule
            dd_alpha_i = (motor_torque_i - resistance_torque_to_motor_i)/self.rotor_inertia
            self.dd_alpha.append(dd_alpha_i)
            d_alpha_i = (self.dd_alpha[-1] + self.dd_alpha[-2])*(self.time[-1] - self.time[-2])/2 + self.d_alpha[-1]
            self.d_alpha.append(d_alpha_i)
            alpha_i = (self.d_alpha[-1] + self.d_alpha[-2])*(self.time[-1] - self.time[-2])/2 + self.alpha[-1]
            self.alpha.append(alpha_i)

            # calculation of the motor current by its characteristic curve
            motor_current_i = self.current_table(d_alpha_i)
            self.motor_current.append(motor_current_i)

            # calculation of the load position, speed and acceleration by gear ratio division
            beta_i = alpha_i/self.gear_ratio
            self.beta.append(beta_i)
            d_beta_i = d_alpha_i/self.gear_ratio
            self.d_beta.append(d_beta_i)
            dd_beta_i = dd_alpha_i/self.gear_ratio
            self.dd_beta.append(dd_beta_i)

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
