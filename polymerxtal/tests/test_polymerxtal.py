"""
Unit and regression test for the polymerxtal package.
"""

# Import package, test suite, and other packages as needed
import polymerxtal
import pytest
import sys

def test_polymerxtal_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "polymerxtal" in sys.modules
