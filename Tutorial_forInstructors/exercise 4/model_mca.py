#!/usr/bin/python
#####################################################################################
# model_mca.py
# Author: Veronica Llorens-Rico
# Version: 1.1
# Date: December 2015
# Description: Exercises on MCA of the Central Carbon Metabolism of E. coli
# Reference: Chassagnole et al, 2002
#####################################################################################


# LOAD MODULES
# ===================================================================================
import numpy as np
from scipy import integrate
import sys
sys.path.insert(0, 'src')  # necessary to work both in PyCharm and terminal
import numdifftools as ndt
# ===================================================================================


# LOAD DATA AND MODEL
# ===================================================================================
# Import differential equations model------------------------------------------------
# Import parameters
from equations_noperturbation import eqs_nopt
from fluxes import fluxcalc
from equations import eqs
from plotSim import plotSimulation
from readData import *
# ===================================================================================


# Exercise 4.1 Calculate the fluxes of each reaction in the model in the equilibrium
# situation.
# ===================================================================================
# Load the parameters and define the equilibrium conditions (initial conditions in
# other exercises.
parameters = params('parameters.txt')
initial_cond = [0.169491,  # cdhap
                0.0988719,  # ce4p
                0.39666,  # cpg2
                2.1185,  # cpg3
                0.00813025,  # cpgp
                0.398943,  # crib5p
                0.11126,  # cribu5p
                0.273191,  # csed7p
                0.138348,  # cxyl5p
                0.597874,  # cf6p
                0.279883,  # cfdp
                0.64958,  # cg1p
                3.46767,  # cg6p
                0.22137,  # cgap
                0.0552527,  # cglex
                2.65434,  # cpep
                0.804284,  # cpg
                2.66951]  # cpyr


fluxes_eq = fluxcalc(initial_cond, parameters)
print fluxes_eq

# ===================================================================================


# Exercise 4.2.
# ===================================================================================
#


# ===================================================================================


# Exercise 4.3.
# ===================================================================================
#


# ===================================================================================