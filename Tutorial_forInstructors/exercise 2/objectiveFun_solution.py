#!/usr/bin/python
########################################################
# objectiveFun.py
# Author: Veronica Llorens-Rico
# Version: 1.2
# Date: December 2015
# Last update: March 2016
# Description: Defines the objective function to
# minimize to estimate model parameters
# Reference: Chassagnole et al, 2002
# Reference: Villaverde et al, 2015
########################################################



def objFun(parameters):
    """
    For every experiment, for every observable on each exp, and for every sample in the experiment, minimize difference
    between the experimental value and the simulated value
    :param parameters: parameters of the simulation
    :return:
    """
    
    # PART 1
    # We can pass experimental data directly within this function, as it is not going to change
    # Experimental data are extracted from Villaverde et al, 2015
    import numpy as np
    from equations import eqs
    from scipy import integrate
    from readData import expData

    dataexp = expData('expValues.txt')

    
    # PART 2
    # first we define the initial conditions, known because at t=0 we have values for everything
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
                    2,  # cglex
                    2.824,  # cpep
                    0.793,  # cpg
                    2.669]  # cpyr

    # we define timespan of the simulation (range for which we have data=303s) and run the simulation
    tspan = np.arange(0, 303, 0.05)
    simResult = integrate.odeint(eqs, initial_cond, tspan, args=(parameters,))

    # we convert the result of the simulation to a dictionary, for indexing purposes
    simResult = np.transpose(simResult)

    simVals = {'cdhap': simResult[0],
               'ce4p': simResult[1],
               'cpg2': simResult[2],
               'cpg3': simResult[3],
               'cpgp': simResult[4],
               'crib5p': simResult[5],
               'cribu5p': simResult[6],
               'csed7p': simResult[7],
               'cxyl5p': simResult[8],
               'cf6p': simResult[9],
               'cfdp': simResult[10],
               'cg1p': simResult[11],
               'cg6p': simResult[12],
               'cgap': simResult[13],
               'cglcex': simResult[14],
               'cpep': simResult[15],
               'cpg': simResult[16],
               'cpyr': simResult[17]}

    
    # PART 3
    # here we define the 'real' objective function which we will minimize
    objFunVal = 0
    for datapoint in dataexp:
        timeval = float(datapoint[0])
        realval = float(datapoint[1])
        species = datapoint[2]
        noiseval = float(datapoint[3])

        # small trick to get the index of the time of the experiment in the simulation data
        match = -1
        timeIndex = 0
        while match == -1:
            if tspan[timeIndex] > (timeval - 0.01) and tspan[timeIndex] < (timeval + 0.01):
                match = 1
            else:
                timeIndex += 1

        simval = simVals[species][timeIndex]

        objFunVal = objFunVal + ((simval - realval) ** 2) / (noiseval ** 2)

    print objFunVal
    return objFunVal
