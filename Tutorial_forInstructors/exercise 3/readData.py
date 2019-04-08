#!/usr/bin/python
########################################################
# readData.py
# Author: Veronica Llorens-Rico
# Version: 1.1
# Date: December 2015
# Description: Defines parameters of the model of the
#   central carbon metabolism of E. coli
# Reference: Chassagnole et al, 2002
########################################################


# Loads parameters from a specified input file (1 parameter per line)
# =======================================================

def params(filetoread):
    """
    Reads parameters from a *.txt file (one parameter per line), in the order pre-specified
    TO DO: we may change this into a dictionary in the future
    :param filetoread:
    :return:
    """
    # read parameters from a file
    par = []
    with open(filetoread) as parameterFile:
        for line in parameterFile:
            a = line.rstrip('\n')
            par.append(float(a))

    return par  # replace by canonical names OR create a dictionary/tuple

# =======================================================



# Loads exp values from a specified input file
# =======================================================

def expData(filetoread):
    """
    Reads experimental values from a *.txt file (one parameter per line)
    TO DO: we may change this into a dictionary in the future
    File columns = 1.Time, 2.Value, 3.Metabolite, 4.standard deviation
    :param filetoread:
    :return:
    """
    # read parameters from a file
    dataLoaded = []
    with open(filetoread) as dataFile:
        for line in dataFile:
            a=line.rstrip('\n').split('\t')
            dataLoaded.append(a)

    return dataLoaded  # replace by canonical names OR create a dictionary/tuple

   # =======================================================