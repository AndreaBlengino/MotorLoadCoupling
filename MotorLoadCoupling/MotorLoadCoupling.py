from .functions import load_data
from .functions import convert_units
from .functions import calculate_interpolation_tables
from .functions import initialize_array
from .functions import fixed_step_time_integration
from .functions import variable_step_time_integration
from .functions import plot_variables


class Coupling:

    def __init__(self):

        self.parameters = {}
        self.variables = {}

    def set_characteristics_curves_files(self, motor_torque_curve_file, motor_current_curve_file, load_torque_curve_file):

        self.parameters['motor_torque_curve_file'] = motor_torque_curve_file
        self.parameters['motor_current_curve_file'] = motor_current_curve_file
        self.parameters['load_torque_curve_file'] = load_torque_curve_file

    def set_inertia(self, rotor_inertia, load_inertia):

        self.parameters['rotor_inertia'] = rotor_inertia
        self.parameters['load_inertia'] = load_inertia

    def set_gearbox(self, load_repetition, gear_ratio, efficiency):

        self.parameters['load_repetition'] = load_repetition
        self.parameters['gear_ratio'] = gear_ratio
        self.parameters['efficiency'] = efficiency

    def set_initial_conditions(self, beta_0, d_beta_0):

        self.parameters['beta_0'] = beta_0
        self.parameters['d_beta_0'] = d_beta_0

    def set_time(self, simulation_time, step_type):

        self.parameters['simulation_time'] = simulation_time
        self.parameters['step_type'] = step_type

    def set_integrator(self, **kargs):

        if self.parameters['step_type'] == 'fixed':
            self.parameters['time_discretization'] = kargs['time_discretization']
        elif self.parameters['step_type'] == 'variable':
            self.parameters['first_step'] = kargs['first_step']
            self.parameters['min_step'] = kargs['min_step']
            self.parameters['max_step'] = kargs['max_step']
            self.parameters['order'] = kargs['order']
        else:
            raise ValueError("step_type must be 'fixed' or 'variable'")

    def start(self):

        load_data(self.parameters)
        convert_units(self.parameters)
        calculate_interpolation_tables(self.parameters)
        initialize_array(self.parameters, self.variables)

        if self.parameters['step_type'] == 'fixed':
            fixed_step_time_integration(self.parameters, self.variables)
        elif self.parameters['step_type'] == 'variable':
            variable_step_time_integration(self.parameters, self.variables)
        else:
            raise ValueError("step_type must be 'fixed' or 'variable'")

    def plotting(self):

        plot_variables(self.variables)
