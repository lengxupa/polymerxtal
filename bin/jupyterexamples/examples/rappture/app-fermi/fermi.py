# ----------------------------------------------------------------------
#  EXAMPLE: Fermi-Dirac function in Python.
#
#  This simple example shows how to use Rappture within a simulator
#  written in Python.
# ======================================================================
#  AUTHOR:  Michael McLennan, Purdue University
#  Copyright (c) 2004-2012  HUBzero Foundation, LLC
#
#  See the file "license.terms" for information on usage and
#  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
# ======================================================================
import Rappture
import sys
from math import *

# open the XML file containing the run parameters
driver = Rappture.library(sys.argv[1])

Tstr = driver.get('input.(temperature).current')
T = Rappture.Units.convert(Tstr, to="K", units="off")

Efstr = driver.get('input.(Ef).current')
Ef = Rappture.Units.convert(Efstr, to="eV", units="off")

print Tstr, T, Efstr, Ef

kT = 8.61734e-5 * T
Emin = Ef - 10*kT
Emax = Ef + 10*kT

E = Emin
dE = 0.005*(Emax-Emin)

# Label the output graph with a title, x-axis label,
# y-axis label, and y-axis units
driver.put('output.curve(f12).about.label','Fermi-Dirac Factor',append=0)
driver.put('output.curve(f12).xaxis.label','Fermi-Dirac Factor',append=0)
driver.put('output.curve(f12).yaxis.label','Energy',append=0)
driver.put('output.curve(f12).yaxis.units','eV',append=0)

while E < Emax:
    f = 1.0/(1.0 + exp((E - Ef)/kT))
    line = "%g %g\n" % (f, E)
    Rappture.Utils.progress(((E-Emin)/(Emax-Emin)*100),"Iterating")
    driver.put('output.curve(f12).component.xy', line, append=1)
    E = E + dE

Rappture.result(driver)
sys.exit()
