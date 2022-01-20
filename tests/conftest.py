from pytest import fixture
from MotorLoadCoupling import Coupling


@fixture(scope = 'function')
def coupling():
    coupling = Coupling()
    return coupling
