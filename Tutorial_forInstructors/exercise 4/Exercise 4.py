
# coding: utf-8

# # Exercise 4. Metabolic control analysis
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
# ## Perform MCA of the Central Carbon Metabolism model
# 
# Metabolic control analysis has a wide range of applications, interestingly in the field of metabolic engineering. By calculating Flux and Concentration control coefficients, we can determine which reactions and enzymes are more sensitive to perturbations.
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
# =====================================================
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', 500)
from scipy import integrate
import sys
sys.path.insert(0, 'src')
import numdifftools as ndt
import warnings
warnings.filterwarnings('ignore')

# LOAD DATA AND MODEL
# =====================================================
# Import differential equations model------------------
# Import parameters
from fluxes import fluxcalc
from readData import *
import roadrunner


# ### Exercise 4.1. Calculate the fluxes of each reaction in the model in the steady state situation. 
# 
# As mentioned in the presentation, fluxes for each reaction are calculated as the difference between the forward and reverse reaction velocities. *J = Vf - Vr*. This is equal to the calculated reaction rates. 

# In[ ]:

# Load the parameters and define the steady state conditions
parameters = params('parameters.txt')
initial_cond = [0.169491,    # cdhap
                0.0988719,   # ce4p
                0.39666,     # cpg2
                2.1185,      # cpg3
                0.00813025,  # cpgp
                0.398943,    # crib5p
                0.11126,     # cribu5p
                0.273191,    # csed7p
                0.138348,    # cxyl5p
                0.597874,    # cf6p
                0.279883,    # cfdp
                0.64958,     # cg1p
                3.46767,     # cg6p
                0.22137,     # cgap
                0.0552527,   # cglex
                2.65434,     # cpep
                0.804284,    # cpg
                2.66951]     # cpyr


# To determine the fluxes, i.e., the reaction rates, under these steady state conditions, we have written a dedicated function **`fluxcalc`**, located in the **`fluxes.py`** file. Take a look at the file.
# 
# The function is very similar to the one describing the system of ODEs of the model. However, this function returns the reaction rates instead of the concentrations of the different metabolites.
# 
# The function requires only two arguments: the initial conditions, and the set of parameters of the model. 

# In[ ]:

# Calculate the fluxe in the steady state
fluxes_eq = fluxcalc(initial_cond, parameters)
print fluxes_eq


# ### Exercise 4.2. Calculate the matrix of Elasticities for all metabolites and all reactions in the model.
# 
# To perform MCA, we will use the Python module **RoadRunner**. This is a module specifically designed to work with SBML systems biology models, and provides tools for loading, simulating, and analyzing these models. 
# 
# This module provides dedicated functions to run Metabolic Control Analysis, facilitating the analysis. As we have previously seen in the tutorial, derivating control coefficients from branched models is not a trivial problem. 
# 
# ------
# 
# First, we will load an SBML version of the model with the **RoadRunner** module.

# In[ ]:

# Import the SBML model of the Central Carbon Metabolism described in Chassagnole et al, 2002
model = roadrunner.RoadRunner("SBML_glycolysis.xml")


# The model is imported as a *RoadRunner* object. After importing the model, we have to ensure that all the metabolite concentrations are in steady state. We can calculate the steady state of the model using the **`steadyState`** function of the **RoadRunner** module. 
# 
# We can then extract the steady state values using the function **`getSteadyStateValues`**

# In[ ]:

# Calculate steady state
model.steadyState()

# Get steady state values
model.getSteadyStateValues()


# Are these values similar to the ones that we used in the previous exercises?
# 
# -----
# 
# Now, we can calculate the enzyme elasticities with respect to all metabolites. Remember that this is calculated by derivating the reaction rate of a specific enzyme with respect to each of the metabolites, and then scaling. 
# 
# This can be done using the **`getScaledElasticityMatrix`** function from the **RoadRunner** module. 

# In[ ]:

# Obtain matrix of elasticities
elasticities=model.getScaledElasticityMatrix()


# Before printing, it is useful to convert the elasticities array into a **Pandas** `DataFrame` object, much more adequate for visualizing here.

# In[ ]:

# Convert the array into a Pandas DataFrame object and print
namesMetabolites=['cdhap','ce4p','cf6p','cfdp','cg1p','cg6p','cgap','cpep','cpg', 'cpg2', 'cpg3', 'cpgp','cpyr','crib5p','cribu5p','csed7p','cxyl5p','cglex']
namesReactionRates=['vALDO','vDAHPS','vDHAP','vE4P','vENO','vEXTER','vG1PAT','vG3PDH','vG6P','vG6PDH','vGAP','vGAPDH','vGLP','vMURSyNTH','vMethSynth','vPDH','vPEP','vPFK','vPG','vPG3','vPGDH','vPGI','vPGK','vPGM','vPGP','vPK','vPPK','vPTS','vR5PI','vRIB5P','vRibu5p','vRu5P','vSED7P','vSynth1','vSynth2','vTA','vTIS','vTKA','vTKB','vTRPSYNTH','vXYL5P','vf6P','vfdP','vpepCxylase','vpg2','vpyr','vrpGluMu','vsersynth']

elasticities=pd.DataFrame(data=elasticities, index=namesReactionRates, columns=namesMetabolites)

print elasticities


# Why is this matrix sparse?

# ### Exercise 4.3. Determine the Flux Control Coefficient for every reaction on the glucose uptake. 
# 
# The glucose uptake is the first reaction of the system, and it is irreversible. Which reactions are the ones increasing the flux of glucose uptake? Which reactions have a negative effect? Consider only the reactions of the network, ignore dilution effects and also the external supply of glucose.
# 
# For this exercise, we also use the **RoadRunner** module. 
# 
# The function to calculate the Flux Control Coefficients is **`getScaledFluxControlCoefficientMatrix`**

# In[ ]:

# Calculate and print Flux Control Coefficients, converting first to a Pandas DataFrame object
FCCmatrix=model.getScaledFluxControlCoefficientMatrix()

FCCdf=pd.DataFrame(data=FCCmatrix, index=namesReactionRates, columns=namesReactionRates)

print FCCdf


# In this matrix, the columns indicate the 'cause', this is the reactions 'modified', and the rows the 'consequence', the effect in the fluxes.
# 
# We can select the row referring to the flux of glucose uptake to better see which are the reactions influencing it

# In[ ]:

# Select the row corresponding to the glucose uptake
glcUptake=FCCdf.loc['vPTS']

print glcUptake


# Ignoring the external glucose supply, 
# 
# Which reactions *increase* the Flux of glucose uptake the most?
# 
# Which ones *decrease* this Flux the most?
# 
# -------
# 
# Test the summation theorem. Do the perturbations affecting the same flux add up to one?

# In[ ]:

# Test the summation theorem by adding all the Flux Control Coefficients of all the enzymes over the same reaction
np.sum(FCCmatrix[[1], ])


# Test the connectivity theorem. Choose any metabolite from the model and multiply the reactions 'connecting' with this metabolite by the elasticities of the enzymes involved. Does this add up to zero?

# In[ ]:

# Test the connectivity theorem


# ### Exercise 4.4. Determine the Concentration Control Coefficient for every reaction on the pyruvate concentration. 
# 
# Imagine that from an engineering point of view, we want to increase the concentration of pyruvate in our system. To do so, we can determine the Concentration Control Coefficients for all reactions on pyruvate, and then modify enzyme concentrations (i.e. overexpress or knock-down) accordingly.
# 
# Calculate the Concentration Control Coefficient matrix using the function **`getScaledConcentrationControlCoefficientMatrix`**.

# In[ ]:

# Calculate and print Concentration Control Coefficients
CCCmatrix=model.getScaledConcentrationControlCoefficientMatrix()

CCCdf=pd.DataFrame(data=CCCmatrix, index=namesMetabolites, columns=namesReactionRates)

print CCCdf


# In this matrix, the columns indicate the 'cause', that is the reactions 'modified', and the rows the 'consequence', the effect in the different metabolite concentrations. Select the row corresponding to the pyruvate concentration.
# 
# Get the row corresponding to the pyruvate concentration to see how it can be affected by the different reaction rates.

# In[ ]:

# Select the row corresponding to the pyruvate concentration
pyruvate=CCCdf.loc['cpyr']

print pyruvate


# Ignoring the external glucose supply,
# 
# What is the reaction with the largest influence in *increasing* pyruvate concentration?
# 
# What is the reaction with the largest influence in *decreasing* pyruvate concentration?

# In[ ]:




# In[ ]:




# In[ ]:



