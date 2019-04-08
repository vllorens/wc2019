
# coding: utf-8

# # Exercise 3. Stability analysis
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
# ## Assess the stability of the steady state described in the previous exercises
# 
# In the previous exercise, we had a system in steady state which we perturbed with the glucose pulse. What happens to the system when it's perturbed? Does it go back to the same steady state? In other words, is this steady state stable? We will assess this in the following exercise.
# 
# ------
# 
# The model was firstly described in [Chassagnole et al, 2002](http://onlinelibrary.wiley.com/doi/10.1002/bit.10288/full). An implementation in Matlab of the model was obtained from [Villaverde et al, 2015](http://gingproc.iim.csic.es/biopredynbench/documentation.html) and recoded using Python.
# 
# ------
# 
# First, we load the modules, data and custom functions that are required for this exercise

# In[ ]:

# LOAD MODULES
# ===================================================================================
import numpy as np
from scipy import integrate
import sys
sys.path.insert(0, 'src')  # necessary to work both in PyCharm and terminal
import numdifftools as ndt
import warnings
warnings.filterwarnings('ignore')

# LOAD DATA AND MODEL
# ===================================================================================
# Import differential equations model------------------------------------------------
# Import parameters
from equations_noperturbation import eqs_nopt
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


# ### Exercise 3.1. Run the model with the estimated parameters, and the glucose perturbation, for a long time.
# 
# Run the model as in the previous exercise between t0=0 and t=1000. Does the model go back to the initial steady state point? Does it diverge? Does it oscillate?

# In[ ]:

# First, we load the parameters and the initial conditions.
parameters = params('parameters.txt')
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
                2,      # cglex
                2.824,  # cpep
                0.793,  # cpg
                2.669]  # cpyr


# Then, we run the model, defining the time span, and plot the results
t = ######<<<>>>######
ds = ######<<<>>>######
plotSimulation(ds, t)


# ### Exercise 3.2. Calculate the Jacobian matrix of the model and solve it for the steady state point. Determine the eigenvalues of this matrix and describe the type of steady state. 
# 
# **Note:** at the steady state, we have the same initial conditions as before, but no glucose perturbation

# In[ ]:

# Modify the initial concentration of external glucose (remember the index is 14)
######<<<>>>######


# As there is no perturbation, we use a modified version of the **`equations.py`** script that does not account for the time dependencies of the NADH, NADP, NADPH, ATP, ADP and AMP. This is the **`equations_noperturbation.py`**, from which the eqs function (**`eqs_nopt`**) has already been imported.
# 
# The function that computes the Jacobian for a system of Ordinary Differential Equations such as the one presented here is **`ndt.Jacobian`**. The only required argument is the equations function.

# In[ ]:

# Calculate the Jacobian for the modified equations
Jacob = ndt.Jacobian(eqs_nopt)


# Now, once the Jacobian matrix has been calculated, we need to pass the steady state concentrations (same as initial conditions) to solve the Jacobian for the steady state point. They need first to be converted to a NumPy array.

# In[ ]:

# Pass steady state concentrations as a NumPy array
j_matrix_solved = Jacob(np.array(initial_cond))


# With the function **`linalg.eig`** from `NumPy` we can get the eigenvalues of a square matrix. This function returns the eigenvalues and the eigenvectors. We are only interested in the first.

# In[ ]:

# Extract the eigenvalues of the Jacobian matrix
eigenvalues = np.linalg.eig(j_matrix_solved)[0]


# Now print the eigenvalues. Can you determine the type of equilibrium? Does any of the eigenvalues have an imaginary part? 
# 
# Does this coincide with your observation from exercise 3.1?

# In[ ]:

print eigenvalues


# ### Exercise 3.3. In the model, we have small dilution effects. How they affect the steady state?
# 
# In the model, we have small dillution effects (D=0.278e-4 s-1). Do they affect how the system responds after a perturbation? 
# 
# The reciprocal absolute values of the real part of the eigenvalues are the time constants of the decline (determining how long takes to the system to go back to the equilibrium). Calculate these values and compare them with the dilution rates.
# 
# ------
# 
# Eigenvalues have a real and an imaginary part. First, we extract the real part using the function **`real`** from NumPy:
# 

# In[ ]:

# Extract the real part of the eigenvalues
eig_real = np.real(eigenvalues)


# The reciprocal of X is calculated as 1/X

# In[ ]:

# Now we find the reciprocals of the real part
reciprocals = abs(1/eig_real)


# Compare these values with the dilution rate. Is dilution affecting how the system behaves?

# In[ ]:

print reciprocals


# In[ ]:



