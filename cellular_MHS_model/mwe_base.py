#!/usr/bin/env python

import os,pdb,re,subprocess,shutil,sys,glob
import json, pickle, time
import itertools, collections, copy
from random import choice
from optparse import OptionParser
from tools import call
from interface import CopasiState
from Truth_Table import GenerateTable, Calculate_State
from ordereddict import OrderedDict
from CopasiRun import OdeModule

def find_changes(Interface):
    """
    Determine the adjusted initial concentration
    of the interface species (ERK, AKT, Raf and
    PTEN) based on their boolean node status
    """
    dict_changes = {}
    #interface_data = dict((key,bool_node_states[key]) for key in sets['thresh'].keys() )
    #Collect interface node data and set changes for the next loop
    #for key, value in interface_data.items():
    #    dict_changes = dict(((v[-1],v[0] + '@' + k),{'value':Peak_Conc[key] * (0.6*value + 0.33*(1 - value))}) for k,v in sets['interface'].items())

    #Loop when all information is passed in a single dictionary
    dict_changes = dict(((value['interface'][-1],value['interface'][0] + '@' + key),{'value':value['PeakConc'] * (0.67*value['status'] + 0.33*(1 - value['status']))}) for key,value in Interface.items())


    changes_next = dict_changes if dict_changes != {} else None
    return changes_next

    

def Combined_Run(args):

    init_state,sets,bool_node_states,reinitialize_ode,output_dirname,copasi_executable,Unconstrained_Nodes = args

    network_prob = 1.0  #Stores the combined network probability
    start_time = time.time()
    #print 'Run started for init state',init_state

    #Change the input states of boolean
    for key, value in zip(Unconstrained_Nodes,init_state):
        bool_node_states[key] = value

    bool_nodes_states_init = copy.deepcopy(bool_node_states)
    #lambda function for getting fingerprints
    fingp_gen = lambda ordered_dict, cstrip = '[]', creplace = (', ',''): str ( [value for value in ordered_dict.values()] ).strip(cstrip).replace(*creplace)

    #from the probability of individual nodes
    #init_state = map(int, sys.argv[1].split())
    #Input_Data = shelve.open('Input_Data.db')
    #Output_Data = shelve.open(os.path.join(output_dirname,'Output_Data' + ('').join(map(str,init_state)) + '.db' ))
    Output_Data = {}

    #sets = Input_Data['sets']
    #Assign correct pathnames to modify in the cps files
    shortcuts = '\n'.join(['path %s = %s'%(key,val) for key,val in sets['paths'].items()])

    #bool_node_states = Input_Data['bool_node_states']

    #copasi_fn = Input_Data['copasi_fn']

    #Create input cps and csv file for each initial condition
    Copasifiles_Modified = copy.deepcopy(sets["Copasifiles"])
    Ode_Return = {}
    Passed_Info = {}
    for key in sets["Copasifiles"]["Original"].keys():
        #Copasifiles_Modified["Original"][key] = os.path.join(output_dirname,sets["Copasifiles"]["Original"][key])
        Copasifiles_Modified["Original"][key] = sets["Copasifiles"]["Original"][key]
        Copasifiles_Modified["Temporary"][key] = os.path.splitext(sets["Copasifiles"]["Temporary"][key])[0] + ('').join(map(str,init_state)) + '.cps'
        Copasifiles_Modified["Output"][key] = os.path.join(output_dirname, os.path.splitext(sets["Copasifiles"]["Output"][key])[0] + ('').join(map(str,init_state)) + '.csv')
        with open(Copasifiles_Modified["Output"][key], 'w') as fp: fp.write(sets['copasi_header'][key]+'\n')
        shutil.copy(Copasifiles_Modified["Original"][key], Copasifiles_Modified["Temporary"][key])
        Active_Nodes, Steady_State_Single, Peak_Conc = OdeModule(Copasifiles_Modified["Temporary"][key],
                sets['thresh'][key],Copasifiles_Modified["Output"][key],copasi_executable,
                changes=sets['changes'][key],shortcuts=shortcuts)
        Ode_Return[key] = {"Active_Nodes" : Active_Nodes, "Steady_State_Single" : Steady_State_Single, "Peak_Conc" : Peak_Conc}
        for key,value in Active_Nodes.items(): 
            bool_node_states[key] = value if value < 2 else bool_node_states[key]
        Passed_Info["init"] = Ode_Return

    #copasi_outfile_orig = str(copasi_outfile)

    #Call to MAPK module
    #Add header line to the copasi output file
    #Call to ode module for initial run


    bool_node_states_ordered = OrderedDict(bool_node_states)
    network_prob_states = {}
    init_state_fingp = str(init_state).strip('[]').replace(', ','')

    boolean_data = dict((key,[val]) for key,val in bool_node_states.items()) #Dictionary for collecting node status
    boolean_state_trn_prob = dict((key,[val]) for key,val in bool_node_states.items()) #Dictionary for collecting node state transition probability
    boolean_data_ordered = OrderedDict(dict((key[0],key[-1]) for key in sets['p53_network'].keys()))
    #Start combined ode and boolean loop
    #------------------------------------------------
    for loop_run in range(sets["STEPS"]):


        #Begin Simulation
        print 'Initial condition', init_state,'Starting loop',loop_run
        #------------------------------------------


        initial_state = fingp_gen(bool_node_states_ordered)
        bool_node_states, state_trn_prob = Calculate_State(sets["p53_network"], bool_node_states)

        for key in bool_node_states.keys(): bool_node_states_ordered[key] = bool_node_states[key]

        final_state = fingp_gen(bool_node_states_ordered)

        for value in state_trn_prob.values(): network_prob *= value
        network_prob_states[init_state_fingp] = network_prob_states.get(init_state_fingp,[]) + [[(initial_state, final_state), network_prob]]
        network_prob = 1.0 #Reinitialize the network probability



        for key, value in bool_node_states.items(): boolean_data[key] = boolean_data.get(key,[]) + [value]
        for key, value in state_trn_prob.items(): boolean_state_trn_prob[key] = boolean_state_trn_prob.get(key,[]) + [value]
        for key in boolean_data.keys(): boolean_data_ordered[key] = boolean_data[key]

        #Get the list of changes
        Interface = {}
        for key, value in sets["interface"].items():
            copy_key = "AKT_PP" if key == "AKT:P:P" else key
            Interface[key] = {"status" : bool_node_states[copy_key], "thresh" : sets["thresh"]["mapk"][copy_key], "PeakConc" : Ode_Return["mapk"]["Peak_Conc"][copy_key], "interface" : value }


        #Calling the mapk module
        changes_from_boolean = find_changes(Interface)

        changes_next = copy.deepcopy(changes_from_boolean)

        #Reset the ode module based on options provided
            
        if reinitialize_ode:
            shutil.copy(Copasifiles_Modified["Original"]["mapk"], Copasifiles_Modified["Temporary"]["mapk"])
            changes_next.update(dict((key,val) for key,val in sets["changes"]["mapk"].items()))
            shutil.copy(Copasifiles_Modified["Original"]["ar"], Copasifiles_Modified["Temporary"]["ar"])
            
        Active_Nodes, Steady_State_Single, Peak_Conc = OdeModule(Copasifiles_Modified["Temporary"]["mapk"],
            sets['thresh']["mapk"],Copasifiles_Modified["Output"]["mapk"],copasi_executable,
            changes=changes_next,shortcuts=shortcuts)
        Ode_Return["mapk"]["Active_Nodes"] = Active_Nodes
        Ode_Return["mapk"]["Steady_State_Single"] = Steady_State_Single
        Ode_Return["mapk"]["Peak_Conc"] = Peak_Conc
        
        if Active_Nodes != {}:
            for key, value in Active_Nodes.items(): 
                bool_node_states[key] = value if value < 2 else bool_node_states[key]

        changes_next = {}
        act_conc = float(Ode_Return["mapk"]["Steady_State_Single"]["AKT_PP"])/10000.0 + float(Ode_Return["mapk"]["Steady_State_Single"]["ErbB2_P"])/10.0

        ## Add the check for p53 concentration
        if bool_node_states["TP53"]:
            p53_conc = 4000.0
        else:
            p53_conc = 0
            
        ar_conc = float(Ode_Return["ar"]["Peak_Conc"]["Androgen_Receptor_Dimer_Cyt"])*(0.67*bool_node_states["AR"] + 0.33*(1 - bool_node_states["AR"]))
            
        changes_next = {
                ("init_species", "cytoplasm@ACT") : {"value" : act_conc},
                ("init_species", "cytoplasm@TP53") : {"value" : p53_conc},
                ("init_species", "cytoplasm@Androgen_Receptor_Dimer_Cyt") : {"value" : ar_conc}
                }
        if reinitialize_ode:
            changes_next.update(dict((key,val) for key,val in sets["changes"]["ar"].items()))

        Active_Nodes, Steady_State_Single, Peak_Conc = OdeModule(Copasifiles_Modified["Temporary"]["ar"],
            sets['thresh']["ar"],Copasifiles_Modified["Output"]["ar"],copasi_executable,
            changes=changes_next,shortcuts=shortcuts)
        Ode_Return["ar"]["Active_Nodes"] = Active_Nodes
        Ode_Return["ar"]["Steady_State_Single"] = Steady_State_Single
        Ode_Return["ar"]["Peak_Conc"] = Peak_Conc
        
        if Active_Nodes != {}:
            for key, value in Active_Nodes.items(): 
                bool_node_states[key] = value if value < 2 else bool_node_states[key]


        Passed_Info["step{0}".format(str(loop_run))] = Ode_Return
        print 'Initial condition', init_state,'Ending loop',loop_run
    #------------------------------------------
    #End of combined loop
    Output_Data['bool_node_states'] = bool_node_states
    Output_Data['boolean_data'] = boolean_data
    Output_Data['boolean_data_ordered'] = boolean_data_ordered
    Output_Data['boolean_state_trn_prob'] = boolean_state_trn_prob
    Output_Data['Passed_Info'] = Passed_Info
    Output_Data['network_prob_states'] = network_prob_states
    #Output_Data.close()
    with open(os.path.join(output_dirname,'Output_Data' + ('').join(map(str,init_state)) + ".pickle"), "wb") as fp: pickle.dump(Output_Data, fp)
    elapsed_time = time.time() - start_time
    print 'Run ended after %f mins for init state'%(elapsed_time/60),init_state

if __name__ == '__main__':

    print 'Main program'
    init_state = [0, 1]
    Constrained_Nodes = ['CDK1','ERK_PP','AKT_PP','Raf_P']
    demo = True
    with open('Model_Input.json','r') as fp: sets = json.load(fp)
    if demo:
        from combined_ode_boolean import demo_changes
        demo_changes(sets)
    Unconstrained_Nodes = list(set(sets["p53_network"].keys()) - set(Constrained_Nodes))
    init_iterator = list(itertools.product([0, 1], repeat=len(Unconstrained_Nodes)))
    bool_node_states = {}
    bool_node_states.update(dict((key,value["initial_state"]) for key,value in sets["p53_network"].items()))
    for key in sets['changes'].keys():
        sets['changes'][key] = dict((tuple(elem[0]),elem[1]) for elem in sets['changes'][key])
    sets['p53_network'] = dict(((key, value['base'],value['initial_state']),value['input_nodes']) for key, value in sets['p53_network'].items())
    reinitialize_ode = False
    output_data = "Test_Data"
    copasi_executable = "Copasi_Executables/64_Bit/CopasiSE"
    args = list(init_state),sets,bool_node_states,reinitialize_ode,output_data,copasi_executable
    Combined_Run(args)
