"""
Unit and regression test for the measure module.
"""

# Import package, test suite, and other packages as needed
import polymerxtal.crystal
import numpy as np
import pytest


# @pytest.mark.slow
def test_calculate_distance():
    """Test that calculate_distance function calculates what we expect."""

    r1 = np.array([0, 0, 0])
    r2 = np.array([0, 1, 0])

    expected_distance = 1

    calculated_distance = polymerxtal.crystal.calculate_distance(r1, r2)

    assert expected_distance == calculated_distance


@pytest.mark.parametrize(
    "r1, r2, r3, expected_angle",
    [
        (np.array([1, 0, 0]), np.array([0, 0, 0]), np.array([0, 1, 0]), 90),
        (
            np.array([np.sqrt(2) / 2, np.sqrt(2) / 2, 0]),
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            45,
        ),
        (np.array([0, 0, -1]), np.array([0, 1, 0]), np.array([1, 0, 0]), 60),
        (
            np.array([np.sqrt(3) / 2, (1 / 2), 0]),
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            30,
        ),
    ],
)
def test_calculate_angle(r1, r2, r3, expected_angle):
    """Test the calculate_angle function"""

    calculated_value = polymerxtal.crystal.calculate_angle(r1, r2, r3, degrees=True)
    assert pytest.approx(expected_angle, abs=1e-2) == calculated_value


@pytest.mark.parametrize(
    "p1, p2, p3, expected_angle",
    [
        (np.array([1, 0, 0]), np.array([0, 0, 0]), np.array([0, 1, 0]), 90),
        (
            np.array([np.sqrt(2) / 2, np.sqrt(2) / 2, 0]),
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            45,
        ),
        (np.array([0, 0, -1]), np.array([0, 1, 0]), np.array([1, 0, 0]), 60),
        (
            np.array([np.sqrt(3) / 2, (1 / 2), 0]),
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            30,
        ),
    ],
)
def test_calculate_angle_many(p1, p2, p3, expected_angle):

    calculated_angle = polymerxtal.crystal.calculate_angle(p1, p2, p3, degrees=True)

    assert expected_angle == pytest.approx(calculated_angle), "%f %f" % (
        calculated_angle,
        expected_angle,
    )


@pytest.mark.parametrize(
    "r1, r2, r3, expected_angle",
    [
        (np.array([1, 0, 0]), np.array([0, 0, 0]), np.array([0, 1, 0]), 90),
        (
            np.array([np.sqrt(2) / 2, np.sqrt(2) / 2, 0]),
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            45,
        ),
        (np.array([0, 0, -1]), np.array([0, 1, 0]), np.array([1, 0, 0]), 60),
        (
            np.array([np.sqrt(3) / 2, (1 / 2), 0]),
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            30,
        ),
    ],
)
def test_calculate_angle_radians(r1, r2, r3, expected_angle):

    calculated_value = polymerxtal.crystal.calculate_angle(r1, r2, r3)
    expected_angle_radians = np.radians(expected_angle)

    assert pytest.approx(expected_angle_radians, abs=1e-2) == calculated_value
