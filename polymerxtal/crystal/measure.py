"""
This module is for functions which perform measurements.
"""

import numpy as np


def calculate_distance(rA, rB):
    """Calculate the distance between two points.

    Parameters
    ----------
    rA, rB : np.ndarray
        The coordinates of each point.

    Returns
    -------
    distance : float
        The distance between the two points.

    Examples
    --------
    >>> r1 = np.array([0, 0, 0])
    >>> r2 = np.array([0, 0.1, 0])
    >>> calculate_distance(r1, r2)
    0.1
    """
    # This function calculates the distance between two points given as numpy arrays.
    if isinstance(rA, np.ndarray) is False or isinstance(rB, np.ndarray) is False:
        raise TypeError("rA and rB must be numpy arrays")

    dist_vec = rA - rB
    distance = np.linalg.norm(dist_vec)

    if distance == 0.0:
        raise Exception("Two atoms are located in the same point in space")

    return distance


def calculate_angle(rA, rB, rC, degrees=False):
    # Calculate the angle between three points. Answer is given in radians by default, but can be given in degrees
    # by setting degrees=True
    AB = rB - rA
    BC = rB - rC
    theta = np.arccos(np.dot(AB, BC) / (np.linalg.norm(AB) * np.linalg.norm(BC)))

    if degrees:
        return np.degrees(theta)
    else:
        return theta


# Conversions between radians and degrees
def RAD2DEG(rad):
    return rad * 57.295779513082325


def DEG2RAD(deg):
    return deg * 0.017453292519943


def calculate_volume(a, b, c, alpha, beta, gamma):
    """Calculate the volume of a unit cell.

    Parameters
    ----------
    a : float
        Unit cell length in Angstroms
    b : float
        Unit cell length in Angstroms
    c : float
        Unit cell length in Angstroms
    alpha : float
        Unit cell angle in degrees
    beta : float
        Unit cell angle in degrees
    gamma : float
        Unit cell angle in degrees

    Returns
    -------
    volume : float
        the volume of the unit cell in cubic Angstrom
    """

    c_x = c * np.cos(DEG2RAD(beta))
    c_y = (
        c
        * (np.cos(DEG2RAD(alpha)) - np.cos(DEG2RAD(beta)) * np.cos(DEG2RAD(gamma)))
        / np.sin(DEG2RAD(gamma))
    )
    c_z = np.sqrt(c * c - c_x * c_x - c_y * c_y)
    return a * b * c_z * np.sin(DEG2RAD(gamma))
