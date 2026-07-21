#!/usr/bin/env python

import os,re,subprocess,shutil,sys,glob
import json, pickle, multiprocessing,platform
import itertools, collections, copy, errno, random
import operator
from random import choice
from optparse import OptionParser
from multiprocessing import Pool

#Local functions
from tools import call
from interface import CopasiState
from Truth_Table import GenerateTable, Calculate_State
from interpret_link import Target_mRNA_List
from ordereddict import OrderedDict
from mwe_base import Combined_Run
from mwe_collector import Collect_Data

def getFromDict(dataDict, mapList):
        return reduce(operator.getitem, mapList, dataDict)

def setInDict(dataDict, mapList, value):
        getFromDict(dataDict, mapList[:-1])[mapList[-1]] = value

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

def MirnaODE(mirt_fn,pat_fn,input_fn,shortcuts,genes):

    """
    This function takes the dictionary of mRNA names
        and changes the corresponding species in ODE model
        by a fixed percentage
    """

    sim = CopasiState(input_fn)
    sim.source_script(shortcuts)

    #Custom expression to remove brackets from an expression
    debracket = lambda x : re.findall('\[([^\]]+)\]',x)
    custom = lambda x : re.findall('^([^\[]+)', x)
    catalog = []
    values = []
    change_list = []
    for elem in sim.get('init_species'):

        species_comp = dict([(key,val[0]) for key,val in [(custom(i.split('=')[1])[0], debracket(i.split('=')[1])) for i in elem.attrib['cn'].split(',')] if val!=[]])
        catalog.append(species_comp)
        values.append(dict({species_comp['Metabolites']:elem.attrib['value']}))

    catalog = [(i['Metabolites'],i['Compartments'],j[i['Metabolites']]) for i,j in zip(catalog,values)]

    for index in [ii for ii,i in enumerate(catalog) if i[0] in genes]:
        sim.change('init_species','%s@%s'%catalog[index][:2],value=float(catalog[index][2])*0.75)
        change_list.append(catalog[index][0])

        sim.write(input_fn,override=True)

    return change_list

def MirnaBool(mirt_fn,pat_fn,p53_network):
    """
    This function takes the dictionary of miRNA names and
    returns the list of target mRNA for each module
    """
    #bool_list = [[key for key,value in i.items() if value == ["bool"]] for i in p53_network.values()]
    #miRNA_Targets_Name_Bool = set([j for i in bool_list for j in i])
    #miRNA_Targets_Bool = dict(zip(miRNA_Targets_Name_Bool, [1]*len(miRNA_Targets_Name_Bool)))

    #ode_list = [[key for key,value in i.items() if value[0] == "ode"] for i in p53_network.values()]
    #miRNA_Target_Name_Ode = set([j for i in ode_list for j in i])
    #miRNA_Targets_Ode = dict(zip(miRNA_Target_Name_Ode,[0.5]*len(miRNA_Target_Name_Ode)))

    miRNA_Targets_Bool, gene_miRNA, gene_all_unpacked = Target_mRNA_List(mirt_fn,pat_fn,p53_network)
    with open('miRNA_List.json','w') as fp: json.dump(gene_miRNA, fp)


    return miRNA_Targets_Bool

def byteify(input):
    """
    This function change json unicode output to ascii
    Source- http://stackoverflow.com/a/13105359
    """

    if isinstance(input, dict):
        return dict((byteify(key),byteify(value)) for key,value in input.iteritems())
    elif isinstance(input, list):
        return [byteify(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input

def Post_Processing(output_data, options):
    """
    Optional post processing of data
    """

    cell_kill_frac, cell_growth_frac,Unconstrained_Nodes = Collect_Data(data_folder=output_data)
    result_string = 'cell_death,cell_growth,cell_senescence\n%f,%f,%f'%(cell_kill_frac,cell_growth_frac,1 - cell_kill_frac - cell_growth_frac)
    with open(os.path.join(output_data, 'cell_fates.csv'),'w') as fp: fp.write(result_string)
    with open(os.path.join(output_data,'ckr_total'),'w') as fp: fp.write('%f'%cell_kill_frac)


    return result_string

if __name__ == '__main__':

    src = os.getcwd()
    curr_dir = os.getenv('PWD')

    #Counters for cell fates
    cell_death = cell_proliferation = cell_senescence = 0

    #Start of command line processing
    #------------------------------------------------

    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.set_defaults(
            mirt_fn = "hsa_MTI.csv",
            Dox=False,
            parallel_run=False,
            no_of_cores=0,
            reinitialize_ode = False,
            output_dirname="Output_Data")

    parser.add_option("-o", "--output-dirname", dest="output_dirname", help="Full or relative(w.r.t. model folder) pathname to the output directory")
    parser.add_option("-M", "--mirtar-file", dest="mirt_fn", help="Name of the mirtarbase data file in csv format")

    parser.add_option("-D", "--Doxorubicin-On", action="store_true", dest="Dox", help="Set whether Doxorubicin is above IC50 or not")
    parser.add_option("-p", "--parallel-run", action="store_true", dest="parallel_run", help="Run program in parallel for different initial conditions, default number of cores is determined by program from available ones")
    parser.add_option("-n", "--no-of-cores", dest="no_of_cores", help="Provide the number of compute cores, by default program determines from available cores")
    parser.add_option("-r", "--reinitialize-ode", action="store_true", dest="reinitialize_ode", help="Specify whether or not to reinitialize the ODE module at every time step of the Boolean(DEFAULT - False)")

    (options, args) = parser.parse_args()

    #End of command line processing

    #Get system architecture and use appropriate copasi file
    if platform.architecture()[0] == '64bit':
        copasi_executable = os.path.join("Copasi_Executables","Linux_64Bit","CopasiSE")
        copasi_executable = 'CopasiSE'
    elif platform.architecture()[0] == '32bit':
        copasi_executable = os.path.join("Copasi_Executables","Linux_32Bit","CopasiSE")
    else:
        Exception("Platform not supported")

    #Get input filename and check existence
    #-----------------------------------------------
    parameter_fn = args[0]

    output_data = options.output_dirname if os.path.isabs(options.output_dirname) else os.path.join(src, options.output_dirname)

    if not os.path.isfile(parameter_fn): raise Exception('[ERROR] invalid parameter file %s'%parameter_fn)
    #Cleanup working directory
    try:
        os.makedirs(output_data)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(output_data): pass
        else: raise
    Input_Data = {}

    with open(parameter_fn,'r') as fp: model_data = byteify(json.load(fp))

    #---custom unpacker
    sets = {}
    for key in model_data: sets[key] = model_data[key]

    #Provide proper extension to provided output filenames
    boolean_outfile = sets["Boolfiles"]
    copasi_outfile = sets["Copasifiles"]["Output"]["mapk"]
    ode_ss_outfile = sets["Copasifiles"]["Ode_SS"]
    input_fn=os.path.join(output_data, sets["Copasifiles"]["Temporary"]["mapk"])


    copasi_outfile_orig = str(copasi_outfile)


    #Add CDK1 and interface species into the list of constrained
    #nodes
    Constrained_Nodes = sets["Constrained_Nodes"]

    #Get and update the initial status of the boolean network
    #------------------------------------------------
    bool_node_states = {}
    bool_node_states.update(dict((key,value["initial_state"]) for key,value in sets["p53_network"].items()))

    #Create lambda for getting the network state
    fingp_gen = lambda ordered_dict, cstrip = '[]', creplace = (', ',''): str ( [value for value in ordered_dict.values()] ).strip(cstrip).replace(*creplace)

    shortcuts = '\n'.join(['path %s = %s'%(key,val) for key,val in sets['paths'].items()])
    #Copy input template copasi file
    copasi_orig_fn = sets["Copasifiles"]["Original"]["mapk"]
    copasi_fn = os.path.join(output_data, os.path.splitext(copasi_orig_fn)[0] + '_adjMirna.cps')
    shutil.copyfile(copasi_orig_fn, copasi_fn)
    #------------------------------------------------

    #Find the list of unconstrained nodes
    Unconstrained_Nodes = sets["Unconstrained_Nodes"]
    print "Prostate Cancer Multiscale Model"
    print "Created by: Dr. Ravi Radhakrishnan, Alok Ghosh"

    #Reorganize the data structure obtained from input file
    #------------------------------------------------
    for key in sets['changes'].keys():
        sets['changes'][key] = dict((tuple(elem[0]),elem[1]) for elem in sets['changes'][key])
    sets['p53_network'] = dict(((key, value['base'],value['initial_state']),value['input_nodes']) for key, value in sets['p53_network'].items())

    Input_Data['sets'] = sets
    Input_Data['bool_node_states'] = bool_node_states
    Input_Data['copasi_fn'] = copasi_fn
    Input_Data['input_fn'] = input_fn
    Input_Data['copasi_outfile'] = copasi_outfile
    Input_Data['boolean_outfile'] = boolean_outfile
    Input_Data['ode_ss_outfile'] = ode_ss_outfile
    Input_Data['Unconstrained_Nodes'] = Unconstrained_Nodes


    jobs = []

    #if options.demo:
    init_iterator = sets["init_iterator"]

    Input_Data['init_iterator'] = init_iterator
    #else:
    #    init_iterator = list(itertools.product([0, 1], repeat=len(Unconstrained_Nodes)))
    #Input_Data.close()
    with open(os.path.join(output_data,'Input_Data.pickle'),"wb") as fp: pickle.dump(Input_Data, fp)


    #Start main loop over all set of initial conditions------------------------------
    #-------------------------------------------------------------------------------

    #---create an argument list for each initial condition
    if options.parallel_run:
        
        joblist = [(list(init_state),sets,bool_node_states,options.reinitialize_ode,output_data,copasi_executable,Unconstrained_Nodes) 
                for init_state in init_iterator]
        #---parallel loop over initial conditions
        if options.no_of_cores:
            pool = Pool(options.no_of_cores)
        else:
            pool = Pool(8)
        results = pool.map(Combined_Run,joblist)
        pool.close()
        pool.join()
    else:

        for init_state in init_iterator:
            
            args = list(init_state),sets,bool_node_states,options.reinitialize_ode,output_data,copasi_executable,Unconstrained_Nodes
            Combined_Run(args)

    result_string = Post_Processing(output_data, options)


    for files in glob.glob("MAPK*.cps"):
        os.remove(files)
    for files in glob.glob("AR_Module_COPASI*.cps"):
        os.remove(files)
    os.chdir(src)
    print "Simulation Results:"
    print result_string
