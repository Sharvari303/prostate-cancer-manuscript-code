#!/usr/bin/env python

import os,pdb,re,subprocess,shutil,sys,glob
import json, shelve
import itertools, collections
from random import choice
from optparse import OptionParser
from tools import call
from interface import CopasiState
from Truth_Table import GenerateTable, Calculate_State
from ordereddict import OrderedDict

def OdeModule(input_fn,thresh,copasi_outfile,copasi_executable,changes=None,shortcuts=None):
    """
    Function for running ODE Network

    This function calls COPASI from within Python using
    subprocess and runs the ODE model for specified duration
    of time. It also returns the active nodes at the end
    of the run based on threshold levels for each species.
    Required arguments -
    input_fn (str) - The name of the input COPASI file
    thresh (dict) - Threshold values of the interface parameters
    copasi_outfile (str) - Name of the output .csv file from Copasi run
    Optional arguments -
    changes (dict) - Contains names and modified values of
    sets of parameters which will be modified in the current
    run
    shortcuts (str) - Contains the full path to navigate to the
    parameter which will be changed

    """

    #Call the CopasiState function from the interface module
    #to parse the input cps file

    sim = CopasiState(input_fn)

    #Define dictionary for steady state values
    Steady_State_Conc = {}
    Peak_Conc = {}

    #Passes the path to the parameters which need to be modified
    #in some manner in this imported xml structure
    #to the source_script function of the interface.py module

    if shortcuts != None: sim.source_script(shortcuts)

    #Implements the changes as specified in the change dictionary
    #argument in the imported xml data structure and overwrites
    #the xml file

    if changes != None:
        for key,val in changes.items(): sim.change(*key,**val)

    #Change the copasi report filename
    change_report = {('report',): {'target': copasi_outfile}}
    for key,val in change_report.items(): sim.change(*key,**val)

    sim.write(input_fn,override=True)

    #---call copasi and run the (modified) ODE model. The .cps
    #file has been setup so that Copasi will update the initial
    #state of species with the values obtained at the end of the run

    output_fn = input_fn
    p = subprocess.Popen(copasi_executable + ' --verbose ' + input_fn + ' -s '+ output_fn,shell=True,
            stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    streams = p.communicate()


    #Get data from Copasi output csv file
    os.system('sed -i "/nan/d" {0}'.format(copasi_outfile))
    with open(copasi_outfile,'r') as fp: time_course_nd = fp.read().splitlines()
    time_course_header = time_course_nd[0].lstrip('#').split(',')

    #Convert the data from the numpy array to a dictionary keyed on headers
    time_course_list = [record.split(',') for record in time_course_nd[1:]]
    time_course_data = dict(zip(time_course_header, [list(t) for t in zip(*time_course_list)]))

    #Get the length of each block from the input cps file (not updated)
    duration = int(sim.get('problem','Duration').attrib['value'])
    stepsize = int(sim.get('problem','StepSize').attrib['value'])
    block_length = duration/stepsize

    #Get a dictionary of continuous concentration from the last block
    continuous_conc = {}
    for key, value in time_course_data.items():
        continuous_conc[key] = map(float, value[len(time_course_list) - int(block_length) : -1])
    #continuous_conc = dict(zip(time_course_data.keys(),[it[len(time_course_data)-int(block_length):-1] for it in time_course_data.values()]))
    Peak_Conc.update(dict((key,max(value)) for key, value in continuous_conc.items()))

    #Get steady state values
    Steady_State_Conc.update(dict((key,value[-1]) for key, value in continuous_conc.items()))
    #Find the discrete concentration and active nodes based
    #on the provided threshold parameters
    discrete_conc,Active_Nodes = {},{}
    for key in thresh:
        #species = sim.get('init_globals','total %s'%key)
        #Distinguish between AR and MAPK cases
        if key == "AR": 
            tempkey = "Androgen_Receptor_Dimer_Cyt"
            if float(continuous_conc[tempkey][-1]) < 0.33*float(Peak_Conc[tempkey]):
                discrete_conc[key] = 0
                Active_Nodes[key] = 0
            elif float(continuous_conc[tempkey][-1]) > 0.6*float(Peak_Conc[tempkey]):
                discrete_conc[key] = 1
                Active_Nodes[key] = 1
            else:
                discrete_conc[key] = 2
                Active_Nodes[key] = 2
        else:
            tempkey = key
            if float(continuous_conc[tempkey][-1]) < 0.33*float(Peak_Conc[tempkey]):
                discrete_conc[key] = 0
                Active_Nodes[key] = 0
            elif float(continuous_conc[tempkey][-1]) > 0.6*float(Peak_Conc[tempkey]):
                discrete_conc[key] = 1
                Active_Nodes[key] = 1
            else:
                discrete_conc[key] = 2
                Active_Nodes[key] = 2

    return Active_Nodes, Steady_State_Conc, Peak_Conc

if __name__ == '__main__':

    print 'Main ODE module'
