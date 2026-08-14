import numpy as np

from icats.functions import EckartAngularDecomposition


def _geometry():
    mass = np.array([16.0, 1.0, 1.0])
    x0 = np.array([[0.0, 0.0, 0.12], [0.9, 0.0, -0.84], [-0.3, 0.82, -0.84]])
    x0 -= np.sum(mass[:, None] * x0, axis=0) / np.sum(mass)
    return x0, mass


def _remove_eckart_leakage(field, x0, mass):
    field = field.copy()
    field -= np.sum(mass[:, None] * field, axis=0) / np.sum(mass)
    inertia = np.zeros((3, 3))
    leakage = np.sum(mass[:, None] * np.cross(x0, field), axis=0)
    for xi, mi in zip(x0, mass):
        inertia += mi * ((xi @ xi) * np.eye(3) - np.outer(xi, xi))
    correction = np.linalg.pinv(inertia) @ leakage
    return field - np.cross(correction, x0)


def _fields():
    x0, mass = _geometry()
    u = _remove_eckart_leakage(
        np.array([[0.02, -0.03, 0.01], [-0.04, 0.05, 0.02], [0.03, 0.01, -0.04]]),
        x0,
        mass,
    )
    du = _remove_eckart_leakage(
        np.array([[-0.03, 0.01, 0.02], [0.04, -0.02, 0.01], [0.02, 0.03, -0.01]]),
        x0,
        mass,
    )
    return x0, x0 + u, du, mass


def test_general_snapshot_decomposition_closes():
    x0, xx, du, mass = _fields()
    omega = np.array([0.17, -0.09, 0.13])
    vv = np.cross(omega, xx) + du

    result = EckartAngularDecomposition(x0, xx, vv, mass, "instantaneous")

    np.testing.assert_allclose(result["omega"], omega, atol=1.0e-12)
    np.testing.assert_allclose(result["internal_velocity"], du, atol=1.0e-12)
    np.testing.assert_allclose(
        result["Iinstant"],
        result["I0"] + result["K1"] + result["K1"].T + result["I2"],
        atol=1.0e-12,
    )
    np.testing.assert_allclose(result["closure"], 0.0, atol=1.0e-12)
    np.testing.assert_allclose(result["eckart_leakage"], 0.0, atol=1.0e-12)


def test_generated_reference_rotation_decomposition_closes():
    x0, xx, du, mass = _fields()
    omega = np.array([0.17, -0.09, 0.13])
    vv = np.cross(omega, x0) + du

    result = EckartAngularDecomposition(x0, xx, vv, mass, "reference")

    np.testing.assert_allclose(result["omega"], omega, atol=1.0e-12)
    np.testing.assert_allclose(result["internal_velocity"], du, atol=1.0e-12)
    np.testing.assert_allclose(result["delta_J2"], 0.0, atol=1.0e-12)
    np.testing.assert_allclose(result["closure"], 0.0, atol=1.0e-12)


def test_decomposition_is_covariant_under_frame_rotation():
    x0, xx, du, mass = _fields()
    omega = np.array([0.17, -0.09, 0.13])
    vv = np.cross(omega, xx) + du
    angle = 0.61
    rotation = np.array(
        [[np.cos(angle), -np.sin(angle), 0.0],
         [np.sin(angle), np.cos(angle), 0.0],
         [0.0, 0.0, 1.0]]
    )

    original = EckartAngularDecomposition(x0, xx, vv, mass, "instantaneous")
    rotated = EckartAngularDecomposition(
        x0 @ rotation.T, xx @ rotation.T, vv @ rotation.T, mass, "instantaneous"
    )

    for key in ("omega", "Jfull", "J0", "delta_J1", "delta_J2", "pi", "closure"):
        np.testing.assert_allclose(rotated[key], rotation @ original[key], atol=1.0e-12)


def test_linear_molecule_uses_nonredundant_rotational_space():
    mass = np.array([14.0, 14.0])
    x0 = np.array([[0.0, 0.0, -0.55], [0.0, 0.0, 0.55]])
    u = np.array([[0.0, 0.0, -0.03], [0.0, 0.0, 0.03]])
    xx = x0 + u
    omega = np.array([0.11, -0.07, 0.0])
    vv = np.cross(omega, xx)

    result = EckartAngularDecomposition(x0, xx, vv, mass, "instantaneous")

    np.testing.assert_allclose(result["omega"], omega, atol=1.0e-12)
    np.testing.assert_allclose(result["internal_velocity"], 0.0, atol=1.0e-12)
    np.testing.assert_allclose(result["closure"], 0.0, atol=1.0e-12)
