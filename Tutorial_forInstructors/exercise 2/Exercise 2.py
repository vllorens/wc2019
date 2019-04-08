
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

# In[ ]:

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
get_ipython().magic(u'matplotlib inline')
matplotlib.rcParams['figure.figsize'] = (15.0, 11.0)
matplotlib.rcParams['font.size'] = 14


# ### Exercise 2.1. A set of initial parameters has been extracted from the literature survey. Run the model with these parameters for a time long enough to  to confirm the system is in a steady state point. Concentrations in the steady state will be our initial conditions. Plot the time courses of the different metabolites. 

# In[ ]:

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

# In[ ]:

# define time limits (use times < 0 as the equations at time 0 are set to introduce a perturbation)
t = np.arange(-50, 0, 0.1)

# integrate
ds = integrate.odeint(eqs, initial_cond, t, args = (initialPar,))

# plot the simulated time-courses of the different metabolites
plotSimulation(ds, t)


# In the previous lines, we used `t<0` as in the model there are some equations that are uncoupled, and have been defined empirically (those for ADP, AMP, ATP, NADH, NADP and NADPH). These are dependent on time. In the steady state, the concentrations of these metabolites are unchanged. In the next exercise, we will apply a perturbation and then we will use `t>0` to incorporate the changes of these metabolites into the model. 

# ### Exercise 2.2. Run the model with these initial conditions from **t0=-10s** to **t=0s**. Run the model again but changing the initial value for external glucose concentration to 2mM (simulating a glucose pulse), from **t0=0s** to **t=350s**. Plot the results of the entire simulation, together with the experimental measurements of 9 metabolites. 
# 
# 
# In the previous exercise, we observed that the initial concentrations were already in steady state as they were not changing over time. Thus we do not need to change the initial parameters. We will run the simulation for 10 seconds and then for 350 seconds more after the perturbation with the addition of glucose to the bioreactor. Experimental data from the timecourse is available for 9 of the metabolites of the model. We will use this data to see how the parameters that we have are able to reproduce the results of the model.
# 
# -----
# 
# First, we run the simulation from `t0=-10s` to `t=0s`, using the previous parameters and initial conditions.

# In[ ]:

# Define time limits for the initial 10s phase
t = ######<<<>>>######

# Integrate
ds = ######<<<>>>######


# Now, run the model again but changing the initial value for external glucose concentration to 2mM (simulating a glucose pulse), from `t0=0s` to `t=350s`. 
# 
# **Note:** The external glucose concentration has index 14 in the `initialPar` array.

# In[ ]:

# Change glucose initial condition for this second simulation
######<<<>>>######

# Define time limits for this second phase of 350s
t2 = ######<<<>>>######

# Integrate
ds2 = ######<<<>>>######


# Now, we have to join the data from the two independent simulations, before and after the perturbation, and plot. To join the two datasets, we use the function **`append`** from NumPy.  
# 
# Before plotting, we load the experimental data using a custom function **`expData`** located in the `readData.py` file. 
# 
# The function **`plotSimulation`** has an optional argument, `experimentalData`, to provide the experimental data read from the file.

# In[ ]:

# join the data from before and after the perturbation and plot
ds = np.append(ds, ds2, axis = 0)
t = np.append(t, t2)

# load the experimental data to be plotted together with the simulations
experimentData = expData('expValues.txt')

# plot the entire simulation with the experimental data
plotSimulation(ds, t, experimentalData = experimentData)


# ### Exercise 2.3. Write an objective function that compares the experimental data measured for 9 metabolites in the model with the simulated data. Use PyMeigo to estimate the set of parameters that best reproduces the experimental data.
# 
# 
# **Meigo** is a global optimization toolbox (see [Egea et al, 2013](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-136)). It has been specifically designed to solve global optimization problems in biology and bioinformatics. It offers different methods that can be used for parameter estimation. One of them (the one that will be used here) is called **enhanced Scatter Search (eSS)**. You can find more about the eSS algorithm [here](http://pubs.acs.org/doi/abs/10.1021/ie801717t).
# 
# Meigo runs in R, but it offers a Python interface called **PyMeigo**. 
# 
# -----
# 
# To run **PyMeigo**, we need to define which parameters need to be optimized. In order to constrain the search, we can define upper and lower bounds for each of the parameters to optimize. The algorithm will not search beyond those limits. We also need to define the *master* function. This is a function that runs the model with some values for the parameters to optimize, evaluates it according to our objective function, and tries to modify those parameters to minimize the value of this objective function.
# 
# -----
# 
# The *master* function has the following parts:
# 
# - **PART 1**: Imports dependent modules and experimental data (explain columns)
# 
# - **PART 2**: Defines the initial conditions (these we know, and will remain unchanged), runs the model with some parameters and gets the results of the simulation. These parameters need to be passed to the master function by PyMeigo, as these are the ones to optimized. 
# 
# - **PART 3**: Objective Function. Evaluates the results of the simulation by comparing to the experimental data. In our case, in the objective function we measure the deviation of the simulated data with the experimental data that we have for 9 of the metabolites. The more experiments, timepoints, and metabolites measured, the more accurate the result will be. 
# 
# The *master* function returns the value of the objective Function, which is what we want to minimize. 
# 
# -----
# 
# For this exercise, we first need to edit the separate script **`objectiveFun.py`**, which contains the master function with the aforementioned three parts. 
# 
# -----
# 
# In **PART 2**, edit the code to run the simulation, considering the adequate timespan. Edit the code to run the simulation with the adequate parameters. 
# 
# -----
# 
# In **PART 3**, for each of the experimental data points, we identify the corresponding species and timepoint in our simulation. Then, the objective function value is updated. Edit the code to update the objective function value correctly, according to the Log-Likelihood cost function.
# 
# <img src="img/daum_equation_1449765695879.png">
# 
# -----
# 
# With these modifications in the **`objectiveFun.py`** code, we can now import this function and run PyMeigo to search for the optimal parameters.

# In[ ]:

# Import the master function
from objectiveFun import objFun  


# Run the optimization with PyMeigo. First, define the upper and lower bounds for the parameters to optimize. Here, the bounds defined are not far from the original values (as they are already pretty optimized), to make the search faster. In real cases, a good choice without prior knowledge would be to multiply the initial values by 0.1 for the lower bound and by 10 for the upper bound.

# In[ ]:

# Define the upper and lower bounds
upperBound = [x*1.5 for x in initialPar]
lowerBound = [x*0.75 for x in initialPar]

# Now, convert the master function to R so that Meigo can run
objfun_for_r=pyfunc2R(objFun)

# Run PyMeigo. It requires only the objective function, and the upper and lower bounds for the parameters to optimize
result_pymeigo = essR(f=objfun_for_r, x_U=upperBound, x_L=lowerBound)


# ### Exercise 2.4. Extract the new set of parameters obtained with either solver and plot the time courses of the 9 metabolites. Add the experimental data points and evaluate the simulations.
# 
# Here, we first extract the set of optimized parameters from the resuts of PyMeigo. Now, we can run the simulation of the perturbation for 350s, as in exercise 2.2, and plot the results. Do these new parameters fit best the data? 

# In[ ]:

# First, we extract the results of PyMeigo as a list
newPar = list(result_pymeigo.xbest)

# Run the simulation of the perturbation, as in exercise 2.2, from t0=0 to t=350s, with the same initial conditions
# (including the glucose pulse)
t = ######<<<>>>######
simulationNew = ######<<<>>>######

# Plot the results of the simulation, with the experimental data
plotSimulation(simulationNew, t, experimentalData=experimentData)


# In the original publication [Chassagnole et al, 2002](http://onlinelibrary.wiley.com/doi/10.1002/bit.10288/epdf), they only plot the datapoints occurring at t<=30s. Compare the results with the ones from the paper. What advantages or disadvantages may have to use only a subset of the data?
# 
# <img src="img/result_paper.png">

# In[ ]:



