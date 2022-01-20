from pytest import mark, raises
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


motor_torque_curve_file = r'data/torque_speed_motor_curve.csv'
motor_current_curve_file = r'data/current_speed_motor_curve.csv'
load_torque_curve_file = r'data/torque_position_load_curve.csv'
rotor_inertia = 4e-4
load_inertia = 1.2
load_repetition = True
gear_ratio = 80
efficiency = 0.7
beta_0 = 0
d_beta_0 = 0
fixed_step_type = 'fixed'
variable_step_type = 'variable'
simulation_time = 50
time_discretization = 0.01
first_step = 1e-6
min_step = 1e-6
max_step = 1e-3
order = 5
tolerance = 1e-6


@mark.set_input
class TestSetInput:

    def test_set_characteristics_curves_files(self, coupling):

        coupling.set_characteristics_curves_files(motor_torque_curve_file = motor_torque_curve_file,
                                                  motor_current_curve_file = motor_current_curve_file,
                                                  load_torque_curve_file = load_torque_curve_file)

        assert coupling.parameters['motor_torque_curve_file'] == motor_torque_curve_file
        assert coupling.parameters['motor_current_curve_file'] == motor_current_curve_file
        assert coupling.parameters['load_torque_curve_file'] == load_torque_curve_file

    def test_set_inertia(self, coupling):

        coupling.set_inertia(rotor_inertia = rotor_inertia, load_inertia = load_inertia)

        assert coupling.parameters['rotor_inertia'] == rotor_inertia
        assert coupling.parameters['load_inertia'] == load_inertia

    def test_set_gearbox(self, coupling):

        coupling.set_gearbox(load_repetition = load_repetition, gear_ratio = gear_ratio, efficiency = efficiency)

        assert coupling.parameters['load_repetition'] == load_repetition
        assert coupling.parameters['gear_ratio'] == gear_ratio
        assert coupling.parameters['efficiency'] == efficiency

    def test_set_initial_conditions(self, coupling):

        coupling.set_initial_conditions(beta_0 = beta_0, d_beta_0 = d_beta_0)

        assert coupling.parameters['beta_0'] == beta_0
        assert coupling.parameters['d_beta_0'] == d_beta_0

    def test_set_time(self, coupling):

        coupling.set_time(step_type = fixed_step_type, simulation_time = simulation_time)

        assert coupling.parameters['step_type'] == fixed_step_type
        assert coupling.parameters['simulation_time'] == simulation_time

    @mark.fixed_step
    def test_set_fixed_integrator(self, coupling):

        coupling.set_time(step_type = fixed_step_type, simulation_time = simulation_time)
        coupling.set_integrator(time_discretization = time_discretization)

        assert coupling.parameters['time_discretization'] == time_discretization

    @mark.variable_step
    def test_set_variable_integrator(self, coupling):

        coupling.set_time(step_type = variable_step_type, simulation_time = simulation_time)
        coupling.set_integrator(first_step = first_step, min_step = min_step, max_step = max_step, order = order)

        assert coupling.parameters['first_step'] == first_step
        assert coupling.parameters['min_step'] == min_step
        assert coupling.parameters['max_step'] == max_step
        assert coupling.parameters['order'] == order

    def test_set_wrong_integrator(self, coupling):

        coupling.set_characteristics_curves_files(motor_torque_curve_file = motor_torque_curve_file,
                                                  motor_current_curve_file = motor_current_curve_file,
                                                  load_torque_curve_file = load_torque_curve_file)
        coupling.set_inertia(rotor_inertia = rotor_inertia, load_inertia = load_inertia)
        coupling.set_gearbox(load_repetition = load_repetition, gear_ratio = gear_ratio, efficiency = efficiency)
        coupling.set_initial_conditions(beta_0 = beta_0, d_beta_0 = d_beta_0)
        coupling.set_time(step_type = 'None', simulation_time = simulation_time)

        with raises(ValueError):
            coupling.set_integrator(time_discretization = time_discretization)

        with raises(ValueError):
            coupling.start()


@mark.analysis
@mark.fixed_step
def test_fixed_step_type(coupling, monkeypatch):

    coupling.set_characteristics_curves_files(motor_torque_curve_file = motor_torque_curve_file,
                                              motor_current_curve_file = motor_current_curve_file,
                                              load_torque_curve_file = load_torque_curve_file)
    coupling.set_inertia(rotor_inertia = rotor_inertia, load_inertia = load_inertia)
    coupling.set_gearbox(load_repetition = load_repetition, gear_ratio = gear_ratio, efficiency = efficiency)
    coupling.set_initial_conditions(beta_0 = beta_0, d_beta_0 = d_beta_0)
    coupling.set_time(step_type = fixed_step_type, simulation_time = simulation_time)
    coupling.set_integrator(time_discretization = time_discretization)

    coupling.start()

    monkeypatch.setattr(plt, 'show', lambda: None)
    coupling.plotting()

    motor_torque_curve = coupling.parameters['motor_torque_curve']
    motor_current_curve = coupling.parameters['motor_current_curve']
    load_torque_curve = coupling.parameters['load_torque_curve']
    motor_torque_table = coupling.parameters['motor_torque_table']
    motor_current_table = coupling.parameters['motor_current_table']
    alpha = np.array(coupling.variables['alpha'])
    d_alpha = np.array(coupling.variables['d_alpha'])
    dd_alpha = np.array(coupling.variables['dd_alpha'])
    beta = np.array(coupling.variables['beta'])
    d_beta = np.array(coupling.variables['d_beta'])
    dd_beta = np.array(coupling.variables['dd_beta'])
    inertia_torque = np.array(coupling.variables['inertia_torque'])
    load_torque = np.array(coupling.variables['load_torque'])
    resistant_torque = np.array(coupling.variables['resistant_torque'])
    resistant_torque_to_motor = np.array(coupling.variables['resistant_torque_to_motor'])
    motor_torque = np.array(coupling.variables['motor_torque'])
    motor_current = np.array(coupling.variables['motor_current'])

    motor_torque_curve_df = pd.read_csv(motor_torque_curve_file)
    motor_torque_curve_df['d_alpha'] = motor_torque_curve_df['d_alpha']*2*np.pi/60
    assert motor_torque_curve_df.equals(motor_torque_curve)

    motor_current_curve_df = pd.read_csv(motor_current_curve_file)
    motor_current_curve_df['d_alpha'] = motor_current_curve_df['d_alpha']*2*np.pi/60
    assert motor_current_curve_df.equals(motor_current_curve)

    load_torque_curve_df = pd.read_csv(load_torque_curve_file)
    load_torque_curve_df['beta'] = load_torque_curve_df['beta']/180*np.pi
    assert load_torque_curve_df.equals(load_torque_curve)

    assert beta[0] == beta_0/180*np.pi
    assert d_beta[0] == d_beta_0*2*np.pi/60
    assert all(alpha/gear_ratio - beta < tolerance)
    assert all(d_alpha/gear_ratio - d_beta < tolerance)
    assert all(dd_alpha/gear_ratio - dd_beta < tolerance)
    assert all(inertia_torque - dd_beta*load_inertia < tolerance)
    assert all(inertia_torque + load_torque - resistant_torque < tolerance)
    assert all(resistant_torque/gear_ratio/efficiency - resistant_torque_to_motor < tolerance)
    assert all(dd_beta - (motor_torque - load_torque/gear_ratio/efficiency)/(gear_ratio*rotor_inertia + load_inertia/gear_ratio/efficiency) < tolerance)
    assert all(motor_torque_table(d_alpha) - motor_torque < tolerance)
    assert all(motor_current_table(d_alpha) - motor_current < tolerance)
    assert all(np.diff(alpha)/time_discretization - d_alpha[1:] < tolerance)
    assert all(np.diff(d_alpha)/time_discretization - dd_alpha[:-1] < tolerance)
    assert all(np.diff(beta)/time_discretization - d_beta[1:] < tolerance)
    assert all(np.diff(d_beta)/time_discretization - dd_beta[:-1] < tolerance)


@mark.analysis
@mark.variable_step
def test_variable_step_type(coupling, monkeypatch):

    coupling.set_characteristics_curves_files(motor_torque_curve_file = motor_torque_curve_file,
                                              motor_current_curve_file = motor_current_curve_file,
                                              load_torque_curve_file = load_torque_curve_file)
    coupling.set_inertia(rotor_inertia = rotor_inertia, load_inertia = load_inertia)
    coupling.set_gearbox(load_repetition = load_repetition, gear_ratio = gear_ratio, efficiency = efficiency)
    coupling.set_initial_conditions(beta_0 = beta_0, d_beta_0 = d_beta_0)
    coupling.set_time(step_type = variable_step_type, simulation_time = simulation_time)
    coupling.set_integrator(first_step = first_step, min_step = min_step, max_step = max_step, order = order)

    coupling.start()

    monkeypatch.setattr(plt, 'show', lambda: None)
    coupling.plotting()

    motor_torque_curve = coupling.parameters['motor_torque_curve']
    motor_current_curve = coupling.parameters['motor_current_curve']
    load_torque_curve = coupling.parameters['load_torque_curve']
    motor_torque_table = coupling.parameters['motor_torque_table']
    motor_current_table = coupling.parameters['motor_current_table']
    alpha = np.array(coupling.variables['alpha'])
    d_alpha = np.array(coupling.variables['d_alpha'])
    dd_alpha = np.array(coupling.variables['dd_alpha'])
    beta = np.array(coupling.variables['beta'])
    d_beta = np.array(coupling.variables['d_beta'])
    dd_beta = np.array(coupling.variables['dd_beta'])
    inertia_torque = np.array(coupling.variables['inertia_torque'])
    load_torque = np.array(coupling.variables['load_torque'])
    resistant_torque = np.array(coupling.variables['resistant_torque'])
    resistant_torque_to_motor = np.array(coupling.variables['resistant_torque_to_motor'])
    motor_torque = np.array(coupling.variables['motor_torque'])
    motor_current = np.array(coupling.variables['motor_current'])

    motor_torque_curve_df = pd.read_csv(motor_torque_curve_file)
    motor_torque_curve_df['d_alpha'] = motor_torque_curve_df['d_alpha']*2*np.pi/60
    assert motor_torque_curve_df.equals(motor_torque_curve)

    motor_current_curve_df = pd.read_csv(motor_current_curve_file)
    motor_current_curve_df['d_alpha'] = motor_current_curve_df['d_alpha']*2*np.pi/60
    assert motor_current_curve_df.equals(motor_current_curve)

    load_torque_curve_df = pd.read_csv(load_torque_curve_file)
    load_torque_curve_df['beta'] = load_torque_curve_df['beta']/180*np.pi
    assert load_torque_curve_df.equals(load_torque_curve)

    assert beta[0] == beta_0/180*np.pi
    assert d_beta[0] == d_beta_0*2*np.pi/60
    assert all(alpha/gear_ratio - beta < tolerance)
    assert all(d_alpha/gear_ratio - d_beta < tolerance)
    assert all(dd_alpha/gear_ratio - dd_beta < tolerance)
    assert all(inertia_torque - dd_beta*load_inertia < tolerance)
    assert all(inertia_torque + load_torque - resistant_torque < tolerance)
    assert all(resistant_torque/gear_ratio/efficiency - resistant_torque_to_motor < tolerance)
    assert all(dd_beta - (motor_torque - load_torque/gear_ratio/efficiency)/(gear_ratio*rotor_inertia + load_inertia/gear_ratio/efficiency) < tolerance)
    assert all(motor_torque_table(d_alpha) - motor_torque < tolerance)
    assert all(motor_current_table(d_alpha) - motor_current < tolerance)
