#!/usr/bin/env python

import os,pdb,re,subprocess,shutil,sys,glob
import json, pickle
import itertools, collections, copy
from random import choice
from optparse import OptionParser
from ordereddict import OrderedDict
from Truth_Table import Calculate_State


def DrugTargets(drug_name,drug_info):
    """
    This function takes data related to drug dosage and
    returns the list of targets in the current network
    for which the initial concentration need to be
    adjusted
    """

    nodes = []

    if drug_name == "Dox":

        #Calculate the atm base level
        KD = drug_info["constants"][0]
        coeff = drug_info["constants"][1]
        atm_base = int(round(3/(((KD/drug_info["dosage"])**coeff) + 1) - 1))

        #Change influenced node status based on threshold
        if drug_info["dosage"] > drug_info["threshold"]:
            nodes = nodes + drug_info["influenced_nodes"]

        return atm_base, nodes

    if drug_name == "vincristine":
        #Check drug threshold value
        cdk1_base = 1
        if drug_info["dosage"] > drug_info["threshold"]:
            nodes = nodes + drug_info["influenced_nodes"]

        return cdk1_base,nodes

def Collect_Data(data_folder='Output_Data',input_filename='Input_Data.pickle'):

    src = os.getcwd()

    if not os.path.isdir(data_folder):
        raise Exception('[ERROR] Output_Data folder does not exists')
    os.chdir(data_folder)
    #Initialize variables
    #------------------------------------------------
    sum_state_prob = 0.0
    count_state_prob = 0
    boolean_SS_data = {} #Dictionary for collecting steady state node status
    boolean_SS_state_prob = {} #Dictionary for collecting steady state node status
    Steady_State_mapk_Combined = {} #Dictionary for collecting steady state concentration from ode module
    Steady_State_ar_Combined = {} #Dictionary for collecting steady state concentration from ode module
    boolean_network_state = {} #Collects the combined state of the network as the collection of numbers
    network_prob_states = {}
    cell_death = cell_proliferation = cell_senescence = 0

    drug_name = 'none'
    drug_list = []

    #Input_Data = shelve.open(input_filename)
    with open(input_filename) as fp: Input_Data = pickle.load(fp)
    sets = Input_Data['sets']
    Chemo_Drug = copy.deepcopy(sets["Chemo_Drug"])

    #print drug_name, src
    for key,val in Chemo_Drug.items():
        if val["status"] == "on": drug_list.append(key)

    if "Dox" in drug_list:
        drug_info = Chemo_Drug["Dox"]
        atm_base,nodes = DrugTargets("Dox",drug_info)
        sets['p53_network']["ATM"]["base"] = atm_base
        for node in nodes:
            bool_node_states[node] = 1
            Constrained_Nodes.append(key)

    if "vincristine" in drug_list:
        drug_info = Chemo_Drug["vincristine"]
        cdk1_base,nodes = DrugTargets("vincristine",drug_info)
        sets['p53_network']["CDK1"]["base"] = cdk1_base


    bool_node_states = Input_Data['bool_node_states']
    copasi_fn = Input_Data['copasi_fn']
    boolean_outfile =  os.path.split(Input_Data['boolean_outfile'])[1]
    ode_ss_outfile_mapk = os.path.split(Input_Data['ode_ss_outfile']["mapk"])[1]
    ode_ss_outfile_ar = os.path.split(Input_Data['ode_ss_outfile']["ar"])[1]
    Unconstrained_Nodes = Input_Data['Unconstrained_Nodes']

    #init_iterator = [map(int, list(re.search('[0-1]+',i).group())) for i in glob.glob('Output*.db')]
    init_iterator = Input_Data["init_iterator"]

    copasi_fn = Input_Data['copasi_fn']
    network_keys  = tuple(OrderedDict(dict((key,value) for key,value in bool_node_states.items())))
    #bool_text_rules = {} #Contains the Boolean state transition rules for each initial state
    #------------------------------------------------

    #Store current boolean node status
    bool_node_states_original = copy.deepcopy(bool_node_states)

    for init_state in init_iterator:

        Init_Temp = zip(Unconstrained_Nodes,list(init_state))
        #Change the initial loop states
        bool_node_states.update(Init_Temp)
        #Generate network state fingerprints
        init_state_fingp = str(init_state).strip('[]').strip('()').replace(', ','')


        #Create input cps and csv file for each initial condition
        input_fn = os.path.splitext(Input_Data['input_fn'])[0] + ('').join(map(str,init_state)) + '.cps'

        copasi_outfile = os.path.splitext(Input_Data['copasi_outfile'])[0] + ('').join(map(str,init_state)) + '.csv'
        copasi_outfile_orig = str(copasi_outfile)

        #Output_Data = shelve.open('Output_Data' + ('').join(map(str,init_state)) + '.db' )
        with open('Output_Data' + ('').join(map(str,init_state)) + '.pickle' ) as fp: Output_Data = pickle.load(fp)

        boolean_data = Output_Data['boolean_data']
        bool_node_states =  Output_Data['bool_node_states']
        boolean_data_ordered = Output_Data['boolean_data_ordered']
        boolean_state_trn_prob = Output_Data['boolean_state_trn_prob']
        Ode_Return = Output_Data['Passed_Info']['step{0}'.format(sets["STEPS"] - 1)]
        Steady_State_mapk_Single = Ode_Return["mapk"]["Steady_State_Single"]
        Steady_State_ar_Single = Ode_Return["ar"]["Steady_State_Single"]
        network_prob_states = Output_Data['network_prob_states']

        #Output_Data.close()


        #Get steady state data for both modules
        #------------------------------------------
        for key, value in boolean_data.items(): boolean_SS_data[key] = boolean_SS_data.get(key,[]) + [value[-1]]
        for key, value in boolean_state_trn_prob.items(): boolean_SS_state_prob[key] = boolean_SS_state_prob.get(key,[]) + [value[-1]]

        for key, value in Steady_State_mapk_Single.items(): Steady_State_mapk_Combined[key] = Steady_State_mapk_Combined.get(key,[]) + [value]
        for key, value in Steady_State_ar_Single.items(): Steady_State_ar_Combined[key] = Steady_State_ar_Combined.get(key,[]) + [value]

        Steady_State_mapk_Combined['Init_FP'] = Steady_State_mapk_Combined.get('Init_FP',[]) + [str(init_state).strip('[]').strip('()').replace(', ','')]
        Steady_State_ar_Combined['Init_FP'] = Steady_State_ar_Combined.get('Init_FP',[]) + [str(init_state).strip('[]').strip('()').replace(', ','')]

        boolean_network_state[network_keys] = boolean_network_state.get(network_keys,[]) + [ [ str ( [value[i] for value in boolean_data_ordered.values()] ).strip('[]').replace(', ','') for i in range ( len ( boolean_data.values()[0] ) ) ] ]
        sum_state_prob  = sum_state_prob + network_prob_states[init_state_fingp][-1][-1]
        count_state_prob = count_state_prob + 1
        #------------------------------------------

        #Get cell death numbers
        #------------------------------------------
        for ind,value in enumerate(boolean_data['CASP9']):
            if not boolean_data["CDK2"][ind] and value:
                cell_death = cell_death + 1
            elif float(Steady_State_ar_Single["PSA"]) > 50.0 or (boolean_data["CDK2"][ind] and not value):

                #When drug is vinchristine and is above threshold
                #value start loop in mitotic arrest
                if drug_name == "vincristine" and nodes != []:
                    bool_node_states["CDK1"] = 1
                    bool_node_states, state_trn_prob = Calculate_State(sets["p53_network"], bool_node_states)

                    if bool_node_states["BCL2"] == 0:
                        cell_death = cell_death + 1

                    if bool_node_states["BCL2"]:
                        cell_proliferation = cell_proliferation + 1
                else:
                        cell_proliferation = cell_proliferation + 1
            else:
                cell_senescence = cell_senescence + 1
        #------------------------------------------

        #Reinitialize the boolean node states
        bool_node_states = copy.deepcopy(bool_node_states_original)


    #------------------------------------------
    #End main loop over all set of initial conditions

    #Get network fingerprint data after loop execution
    #------------------------------------------
    network_fingerprints = dict(("nodes",key) for key in boolean_network_state.keys())
    network_fingerprints.update(dict(("status",value) for value in boolean_network_state.values()))
    #------------------------------------------

    cell_kill_frac = cell_death/((cell_death + cell_proliferation + cell_senescence)*1.0)
    cell_growth_frac = cell_proliferation/((cell_death + cell_proliferation + cell_senescence)*1.0)
    cell_kill_prob = sum_state_prob/count_state_prob

    #Write output files
    #------------------------------------------
    #Check if directory exists
    with open(boolean_outfile,'w') as fp: json.dump(boolean_SS_data, fp)
    with open('State_Transitions_Probability.json','w') as fp: json.dump(boolean_state_trn_prob, fp)
    with open('State_Transitions_Step_Probability.json','w') as fp: json.dump(network_prob_states, fp)
    with open(ode_ss_outfile_mapk,'w') as fp: json.dump(Steady_State_mapk_Combined, fp)
    with open(ode_ss_outfile_ar,'w') as fp: json.dump(Steady_State_ar_Combined, fp)
    with open('State_Transitions.json','w') as fp: json.dump(network_fingerprints, fp)
    with open('Cell_Fate.json','w') as fp: json.dump({'cell_death': [cell_death, cell_kill_frac],'cell_proliferation': [cell_proliferation, cell_growth_frac],'cell_senescence': [cell_senescence, 1 - cell_kill_frac - cell_growth_frac]}, fp)
    #------------------------------------------

    #Calculate cell kill rate
    os.chdir(src)

    return cell_kill_frac, cell_growth_frac, Unconstrained_Nodes

if __name__ == '__main__':

    print 'Main Program'
