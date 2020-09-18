"""
Testing for molecule module
"""

import polymerxtal
import pytest
import sys

import numpy as np
import os

@pytest.fixture(scope="module")
def methane_molecule():
    symbols = np.array(['C', 'H', 'H', 'H', 'H'])
    coordinates = np.array([[1,1,1], [2.4,1,1], [-0.4, 1, 1], [1, 1, 2.4], [1, 1, -0.4]])
    return symbols, coordinates

def test_build_bond_list(methane_molecule):
    symbols, coordinates = methane_molecule

    bonds = polymerxtal.build_bond_list(coordinates)

    assert len(bonds) == 4

    for atoms, bonds in bonds.items():
        assert bonds == 1.4

def test_molecular_mass(methane_molecule):
    symbols, coordinates = methane_molecule
    
    calculated_mass = polymerxtal.calculate_molecular_mass(symbols)

    actual_mass = polymerxtal.atom_data.atomic_weights['C'] + polymerxtal.atom_data.atomic_weights['H'] +\
         polymerxtal.atom_data.atomic_weights['H'] + polymerxtal.atom_data.atomic_weights['H'] + polymerxtal.atom_data.atomic_weights['H']

    assert actual_mass == calculated_mass

def test_build_bond_list_failure():
    coordinates = np.array([])
    
    with pytest.raises(ValueError):
        polymerxtal.build_bond_list(coordinates)

def test_center_of_mass(methane_molecule):
    symbols, coordinates = methane_molecule

    center_of_mass = polymerxtal.calculate_center_of_mass(symbols, coordinates)

    expected_center = np.array([1,1,1])
    
    assert np.array_equal(center_of_mass, expected_center)
