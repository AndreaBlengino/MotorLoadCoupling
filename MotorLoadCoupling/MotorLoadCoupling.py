from .functions import load_data
from .functions import convert_units
from .functions import calculate_interpolation_tables
from .functions import initialize_arrays
from .functions import fixed_step_time_integration
from .functions import variable_step_time_integration
from .functions import plot_variables


class Coupling:
    """
    Time integration of a motor-gearbox-load system.

    Attributes:
        parameters (dict): Dictionary where input parameters of the system are stored.
        variables (dict): Dictionary where time-dependent variables calculated during time integration are stored.

    Methods:
        set_characteristics_curves_files: Sets motor and load's characteristics curves file path.
        set_inertia: Sets rotor and load's moments of inertia.
        set_gearbox: Sets gearbox's parameters.
        set_initial_conditions: Sets load's initial conditions.
        set_time: Sets simulation time parameters.
        set_integrator: Sets integrator's parameters.
        start: Runs simulation.
        plotting: Plots calculated variables.
    """

    def __init__(self):
        """
        Initializes parameters and variables dictionaries.
        The parameters dictionary is used to store input parameters of the system.
        The variables dictionary is used to store time-dependent variables calculated during time integration.
        """

        self.parameters = {}
        self.variables = {}

    def set_characteristics_curves_files(self, motor_torque_curve_file, motor_current_curve_file, load_torque_curve_file):
        """
        Sets motor and load's characteristics curves file path.
        Saves the file path as key-value pairs in the parameters dictionary.
        Files must be in .csv format.
        Three characteristics curves files are required:

        - motor_torque_curve_file reports the driving torque delivered by the motor as a function of its rotation speed.
          The file is structured in two columns:
          * d_alpha which reports the motor rotation speed, expressed in rpm
          * torque which reports the motor torque, expressed in Nm

        - motor_current_curve_file reports the electric current absorbed by the motor as a function of its rotation
          speed.
          The file is structured in two columns:
          * d_alpha which reports the motor rotation speed, expressed in rpm
          * current which reports the motor torque, expressed in A

        - load_torque_curve_file reports the resistant torque exerted by the load as a function of its position.
          The file is structured in two columns:
          * beta which reports the position of the load, expressed in °
          * torque which shows the resistant torque of the load, expressed in Nm

        Args:
            motor_torque_curve_file (str): The file path to the torque-speed motor characteristic curve file.
            motor_current_curve_file (str): The file path to the current-speed motor characteristic curve file.
            load_torque_curve_file (str): The file path to the torque-position load characteristic curve file.

        Returns:
            None.
        """

        self.parameters['motor_torque_curve_file'] = motor_torque_curve_file
        self.parameters['motor_current_curve_file'] = motor_current_curve_file
        self.parameters['load_torque_curve_file'] = load_torque_curve_file

    def set_inertia(self, rotor_inertia, load_inertia):
        """
        Sets rotor and load's moments of inertia.
        Saves them as key-value pairs in the parameters dictionary.

        Args:
            rotor_inertia (float): The rotor's moment of inertia expressed in kgm2.
            load_inertia (float): The load's  moment of inertia expressed in kgm2.

        Returns:
            None.
        """

        self.parameters['rotor_inertia'] = rotor_inertia
        self.parameters['load_inertia'] = load_inertia

    def set_gearbox(self, load_repetition, gear_ratio, efficiency):
        """
        Sets gearbox's parameters.
        Saves them as key-value pairs in the parameters dictionary.

        Args:
            load_repetition (bool): If, during the simulation, the value of the load position beta exceeds the range of
            values​expressed in the load_torque_curve_file file it is necessary to specify what the resistant torque
            value should be.
            If load_repetition is True, then the resistant torque expressed in load_torque_curve_file is repeated
            cyclically as the beta changes; if, on the other hand, it is False, then the resistant torque does not
            repeat itself and the resistant torque value is calculated through a linear extrapolation.
            gear_ratio (float): The total transmission ratio, expressed as the ratio of the rotation speed of the motor
            and the one of the load.
            efficiency (float): The overall efficiency of the gearbox.

        Returns:
            None.
        """

        self.parameters['load_repetition'] = load_repetition
        self.parameters['gear_ratio'] = gear_ratio
        self.parameters['efficiency'] = efficiency

    def set_initial_conditions(self, beta_0, d_beta_0):
        """
        Sets load's initial conditions.
        Saves them as key-value pairs in the parameters dictionary.

        Args:
            beta_0 (float): The initial position of the load, expressed in °.
            d_beta_0 (float): The initial speed of the load, expressed in rpm.

        Returns:
            None.
        """

        self.parameters['beta_0'] = beta_0
        self.parameters['d_beta_0'] = d_beta_0

    def set_time(self, simulation_time, step_type):
        """
        Sets simulation time parameters.
        Saves them as key-value pairs in the parameters dictionary.

        Args:
            simulation_time (float): The total duration of the simulation, expressed in s.
            step_type (str): Indicates whether to use a fixed or variable time step. Can be only 'fixed' or 'variable'.

        Returns:
            None.
        """

        self.parameters['simulation_time'] = simulation_time
        self.parameters['step_type'] = step_type

    def set_integrator(self, **kargs):
        """
        Sets integrator's parameters.
        Saves them as key-value pairs in the parameters dictionary.

        Args:
            time_discretization (float): The time step, expressed in s. Used only if step_type is set to 'fixed'.
            first_step (float): The size of the first time step in s. Used only if step_type is set to 'variable'.
            min_step (float):  The minimum time step in s. Used only if step_type is set to 'variable'.
            max_step (float): The maximum time step in s. Used only if step_type is set to 'variable'.
            order (int): The maximum order used by the integrator, <= 5. Used only if step_type is set to 'variable'.

        Returns:
            None.

        Raises:
            ValueError if step_type is neither 'fixed' nor 'variable'.
        """

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
        """
        Runs simulation:
        - loads characteristics curves from file paths provided through set_characteristics_curves_files
        - converts units
        - calculates interpolation tables for the characteristics curves
        - initializes arrays
        - runs time integration

        Positions, velocities, accelerations, torques and current are saved in variables dictionary.

        Returns:
            None.
        Raises:
            ValueError if step_type is neither 'fixed' nor 'variable'.
        """

        load_data(self.parameters)
        convert_units(self.parameters)
        calculate_interpolation_tables(self.parameters)
        initialize_arrays(self.parameters, self.variables)

        if self.parameters['step_type'] == 'fixed':
            fixed_step_time_integration(self.parameters, self.variables)
        elif self.parameters['step_type'] == 'variable':
            variable_step_time_integration(self.parameters, self.variables)
        else:
            raise ValueError("step_type must be 'fixed' or 'variable'")

    def plotting(self):
        """
        Plots calculated variables against time in a 9x9 grid.

        Returns:
            None.
        """

        plot_variables(self.variables)
