#!/usr/bin/python
########################################################
# plotSim.py
# Author: Veronica Llorens-Rico
# Version: 1.1
# Date: December 2015
# Description: Plots results of the simulation
# Reference: Chassagnole et al, 2002
########################################################


def plotSimulation(simResult, time, experimentalData=None):
    """
    Reformats the result of a simulation and generates a plot of
    the time courses of some of the metabolites
    :param simResult: result of the odeint function
    :param time: time of the integration
    :param experimentalData: boolean stating whether to import and plot the experimental measurements
    :return: plots a time course of some of the metabolites with the simulation data
    """
    import matplotlib.pyplot as plt
    import numpy as np

    simResult = np.transpose(simResult)

    cdhap = simResult[0]
    ce4p = simResult[1]
    cpg2 = simResult[2]
    cpg3 = simResult[3]
    cpgp = simResult[4]
    crib5p = simResult[5]
    cribu5p = simResult[6]
    csed7p = simResult[7]
    cxyl5p = simResult[8]
    cf6p = simResult[9]
    cfdp = simResult[10]
    cg1p = simResult[11]
    cg6p = simResult[12]
    cgap = simResult[13]
    cglcex = simResult[14]
    cpep = simResult[15]
    cpg = simResult[16]
    cpyr = simResult[17]


    # ============================================================
    # Generate plot
    plt.figure()
    p1 = plt.subplot(331)
    plt.plot(time, cglcex, '0.75')
    plt.xlabel('time (s)')
    plt.ylabel('concentrations mM')
    plt.legend(('GLCext',), loc=1, prop={'size' : 10})
    p2 = plt.subplot(332, sharey=p1)
    plt.plot(time, cfdp, '-g')
    plt.xlabel('time (s)')
    plt.ylabel('concentrations mM')
    plt.legend(('FDP',), loc=1, prop={'size' : 10})
    p3 = plt.subplot(333, sharey=p1)
    plt.plot(time, cg1p, '-b')
    plt.xlabel('time (s)')
    plt.ylabel('concentrations mM')
    plt.legend(('G1P',), loc=1, prop={'size' : 10})
    p4 = plt.subplot(334)
    plt.plot(time, cg6p, '-r')
    plt.xlabel('time (s)')
    plt.ylabel('concentrations mM')
    plt.legend(('G6P',), loc=1, prop={'size' : 10})
    p5 = plt.subplot(335, sharey=p4)
    plt.plot(time, cpep, '-c')
    plt.xlabel('time (s)')
    plt.ylabel('concentrations mM')
    plt.legend(('PEP',), loc=1, prop={'size' : 10})
    p6 = plt.subplot(336, sharey=p4)
    plt.plot(time, cpyr, '-m')
    plt.xlabel('time (s)')
    plt.ylabel('concentrations mM')
    plt.legend(('PYR',), loc=1, prop={'size' : 10})
    p7 = plt.subplot(337)
    plt.plot(time, cf6p, '-y')
    plt.xlabel('time (s)')
    plt.ylabel('concentrations mM')
    plt.legend(('F6P',), loc=1, prop={'size' : 10})
    p8 = plt.subplot(338, sharey=p7)
    plt.plot(time, cgap, '-k')
    plt.xlabel('time (s)')
    plt.ylabel('concentrations mM')
    plt.legend(('GAP',), loc=1, prop={'size' : 10})
    p9 = plt.subplot(339, sharey=p7)
    plt.plot(time, cpg, '0.25')
    plt.xlabel('time (s)')
    plt.ylabel('concentrations mM')
    plt.legend(('6PG',), loc=1, prop={'size' : 10})


    if experimentalData is not None:
        for datapoint in experimentalData:
            timeval = float(datapoint[0])
            realval = float(datapoint[1])
            species = datapoint[2]
            if species == 'cpep':
                p5.plot(timeval, realval, 'co')
            elif species == 'cglcex':
                p1.plot(timeval, realval, 'ko')
            elif species == 'cg6p':
                p4.plot(timeval, realval, 'ro')
            elif species == 'cpyr':
                p6.plot(timeval, realval, 'mo')
            elif species == 'cf6p':
                p7.plot(timeval, realval, 'yo')
            elif species == 'cg1p':
                p3.plot(timeval, realval, 'bo')
            elif species == 'cpg':
                p9.plot(timeval, realval, 'ko')
            elif species == 'cfdp':
                p2.plot(timeval, realval, 'go')
            elif species == 'cgap':
                p8.plot(timeval, realval, 'ko')

    plt.tight_layout()
    plt.show()
