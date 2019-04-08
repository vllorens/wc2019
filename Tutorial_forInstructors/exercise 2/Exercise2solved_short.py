
# coding: utf-8

# # Exercise 2. Parameter estimation in ODE models
# 
# Tutorial: Pathway modeling using Ordinary Differential Equations
# 
# Instructor: [Veronica Llorens-Rico](veronica.llorens@crg.eu)
# 
# Whole-cell Model Course; 5th April 2016
# 
# Centre for Genomic Regulation, Barcelona
# 
# 
# ## Write the objective function for this model in Python, and use a solver to estimate the parameters of the model
# 
# In this exercise, we will start with some set of parameters derived from the literature, as in the original publication, and we will optimize these parameters to fit the experimental data as good as possible. Due to time constraints, we will only perform the optimization of a small set of parameters, whilst in the original publication all parameters are optimized.
# 
# The model was firstly described in [Chassagnole et al, 2002](http://onlinelibrary.wiley.com/doi/10.1002/bit.10288/full). An implementation in Matlab of the model was obtained from [Villaverde et al, 2015](http://gingproc.iim.csic.es/biopredynbench/documentation.html) and recoded using Python.
# 
# ------
# 
# First, we load the modules, data and custom functions that are required for this exercise

# In[1]:

# LOAD MODULES
# =====================================================
import numpy as np
from scipy import integrate
from pymeigo import *
import sys
sys.path.insert(0, 'src')


# LOAD DATA AND MODEL
# =====================================================
# Import differential equations model------------------
# Import parameters
from equations import eqs
from plotSim import plotSimulation
from readData import *


# INLINE PLOTS
# =====================================================
import pylab
import matplotlib
import matplotlib.pyplot as plt
#get_ipython().magic(u'matplotlib inline')
#matplotlib.rcParams['figure.figsize'] = (15.0, 11.0)
#matplotlib.rcParams['font.size'] = 14


# ### Exercise 2.1. A set of initial parameters has been extracted from the literature survey. Run the model with these parameters for a time long enough to confirm the system is in a steady state point. Concentrations in the steady state will be our initial conditions. Plot the time courses of the different metabolites. 

# In[2]:

# Load the initial parameters and the initial conditions
initialPar = params('parameters.txt')
initial_cond = [0.185,  # cdhap
                0.103,  # ce4p
                0.422,  # cpg2
                2.254,  # cpg3
                0.008,  # cpgp
                0.393,  # crib5p
                0.108,  # cribu5p
                0.246,  # csed7p
                0.137,  # cxyl5p
                0.570,  # cf6p
                0.334,  # cfdp
                0.616,  # cg1p
                3.307,  # cg6p
                0.242,  # cgap
                0.043,  # cglex
                2.824,  # cpep
                0.793,  # cpg
                2.669]  # cpyr


# **Note:** In order to simulate a system of differential equations over time, we need to define the time limits using the function **`arange`** from NumPy. This function requires the following arguments:
# 
# - Initial time
# 
# - Final time
# 
# - Timestep used to print the results
# 
# To run the simulation, we use the function **`integrate.odeint`**. This function requires the following arguments:
# 
# - Set of Ordinary Differential Equations, as those from the Exercise 1.
# 
# - Set of initial conditions
# 
# - Time limits
# 
# - Arguments (parameters) required by the set of ODEs. In this case, the initial parameters extracted from the literature.
# 
# After the simulation has been executed, we can plot the results by using a custom defined function, **`plotSimulation`**. This function requires the output of the simulation, as well as the integration time.

# In[3]:

# define time limits (use times < 0 as the equations at time 0 are set to introduce a perturbation)
t = np.arange(-50, 0, 0.1)

# integrate
ds = integrate.odeint(eqs, initial_cond, t, args = (initialPar,))

# plot the simulated time-courses of the different metabolites
plotSimulation(ds, t)



