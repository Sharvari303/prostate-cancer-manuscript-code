#!/usr/bin/env python

import os, subprocess, platform,re, itertools
import json, collections, shelve, glob
import xml.etree.ElementTree as ET
import numpy as np

def makedirectory(dirname="New_Folder"):
    """
    Creates a new directory if it
    does not exists
    """

    status = 0
    try:
        os.makedirs(dirname)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            return status
        else:
            raise
    finally:
        status = 1
        return status

def GetCpsvalues(input_cps, ModelParameterGroup="Species"):

    """
    This function returns values from specific Model
    Parameter Groups in a copasi file in the format
    'name', 'value', 'compartment' etc
    """
    #Get the initial concentrations from the cps file
    with open(input_cps, 'r') as fp: xmlstring = fp.read()
    xmlstring = re.sub(' xmlns="[^"]+"','',xmlstring,count=1)
    root = ET.fromstring(xmlstring)

    param_list = []
    change_list = {}
    if ModelParameterGroup == "Species":
        parametergroup_string = 'String=Initial Species Values'
    else:
        parametergroup_string = 'String=Initial Global Quantities'
    
    Model_ParameterGroups = root.findall(".//ModelParameterGroup")
    
    for parametergroup in Model_ParameterGroups:
        if parametergroup.get('cn') == parametergroup_string:
            param_elements = parametergroup.findall("ModelParameter")
            break
        
    else:
        raise ValueError('No such attribute Initial Species Values')

    for elem in param_elements:
        
        if ModelParameterGroup == "Species":
            param_name = re.search('Metabolites\[(.+)\]', elem.get('cn').split(',')[3]).group(1)
            param_compartment = re.search('Compartments\[(.+)\]', elem.get('cn').split(',')[2]).group(1)
            param_value = elem.get('value')
            param_list.append([param_name, param_value, param_compartment])
        else:
            param_name = re.search('Values\[(.+)\]', elem.get('cn').split(',')[2]).group(1)
            param_value = elem.get('value')
            param_list.append([param_name, param_value])

    
    return param_list

def Change_CPS(root, change_set,  compartments=[], copasi_outfile=''):
    """
    This function takes copasi xml root object as
    an argument and changes a single parameter
    in the root 
    """
    
    model_names = root.findall("Model")
    model_name = model_names[0].get('name')


    
    species_string = "CN=Root,Model={0},Vector=Compartments[{1}],Vector=Metabolites[{2}]"
    globals_string = "CN=Root,Model={0},Vector=Values[{1}]"
    change_list = []
    
    if compartments != []:
        parametergroup_string = 'String=Initial Species Values'
        search_string = species_string.format(model_name,compartments[1],change_set[0])
    else:
        parametergroup_string = 'String=Initial Global Quantities'
        search_string = globals_string.format(model_name, change_set[0])

    #print search_string


    Model_ParameterGroups = root.findall(".//ModelParameterGroup")
    for parametergroup in Model_ParameterGroups:
        if parametergroup.attrib.get('cn') == parametergroup_string:
            ModelParameter = parametergroup.findall("ModelParameter")
            for Parameter in ModelParameter:
                if Parameter.attrib.get('cn') == search_string:
                    current_value = Parameter.get('value')
                    change_list = change_list + [change_set[0],current_value]
                    change_to_value = str(change_set[1])
                    Parameter.set('value',change_to_value)
                    change_list.append(Parameter.get('value'))
                    break
            else:
                raise ValueError('No such parameter {0}'.format(change_set[0]))
            break
    else:
        raise ValueError('No such attribute Initial Species Values')

    if copasi_outfile != '':
        listoftasks = root.findall(".//Task")
        for tasks in listoftasks:
            if tasks.get('name') == "Time-Course":
                reports = tasks.find("Report")
                reports.set('target',copasi_outfile)
                break
        else:
            raise ValueError('Report Time-Course not found')

    return change_list




def OdeRun(input_cps, joblimit=4):
    """
    Wrapper function to call Copasi executable and run for
    a specified set of conditions
    """

    #Get system architecture and use appropriate copasi file
    out = ''
    err = ''
    if platform.architecture()[0] == '64bit':
        Copasi_Command = './Copasi_Executables/64_Bit/CopasiSE'
    if platform.architecture()[0] == '32bit':
        Copasi_Command = './Copasi_Executables/32_Bit/CopasiSE'

    #Call Copasi executable using the name of the input cps
    output_cps = input_cps
    p = subprocess.Popen(Copasi_Command + ' --verbose ' + input_cps,shell=True,
            stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    for lines in iter(p.stdout.readline,''): out += lines.rstrip()
    for lines in iter(p.stderr.readline,''): err += lines.rstrip()
    streams = p.communicate()

    return out, err

def Scan_Setup(input_cps, scanfile="Parameter_Scan.json", output_cps='', combined=False):
    """
    This sets up a parameter scan using values set in
    the scan json file
    """
    with open(scanfile,'r') as fp: Scan_Data = json.load(fp)
    change_list = []
    copasi_outfile='Time_Course.csv'
    output_dirname='ScanDir_Species'

    with open(input_cps, 'r') as fp: xmlstring = fp.read()
    xmlstring = re.sub(' xmlns="[^"]+"','',xmlstring,count=1)
    root = ET.fromstring(xmlstring)


    if not combined:
        for key, value in Scan_Data["species"].items():
            dirname = output_dirname + key
            status = makedirectory(dirname=dirname)
            for i,val in enumerate(value[0]):
                value_list = [key,val]
                output_cps = os.path.splitext(input_cps)[0] + str(i) + '.cps'
                output_cps = os.path.join(dirname, output_cps)
                copasi_outfile_num = os.path.splitext(copasi_outfile)[0] + str(i) + '.csv'
                compartment_list = [key,value[1]]
                change_list += Change_CPS(root, value_list,compartments=compartment_list,copasi_outfile=copasi_outfile_num)
                ET.ElementTree(root).write(output_cps, encoding='utf-8')
                root = ET.fromstring(xmlstring)
    else:
        keylist = Scan_Data['species'].keys()
        dirname = output_dirname + '_'.join(keylist)
        status = makedirectory(dirname=dirname)
        combined_value_list = list(itertools.product(*[Scan_Data['species'][key][0] for key in keylist]))
        combined_pair_list = [[[keylist[i], elem[i]] for i in range(len(keylist))] for elem in combined_value_list]
        for i,val in enumerate(combined_pair_list):
            output_cps = os.path.splitext(input_cps)[0] + str(i) + '.cps'
            output_cps = os.path.join(dirname, output_cps)
            copasi_outfile_num = os.path.splitext(copasi_outfile)[0] + str(i) + '.csv'

            for j in range(len(val)):
                value_list = val[j]
                compartment_list = [val[j][0],Scan_Data['species'][val[j][0]][1]]
                change_list += Change_CPS(root, value_list,compartments=compartment_list,copasi_outfile=copasi_outfile_num)
            ET.ElementTree(root).write(output_cps, encoding='utf-8')
            root = ET.fromstring(xmlstring)

    for key, value in Scan_Data["param"].items():
        dirname = output_dirname + key
        status = makedirectory(dirname=dirname)
        for i, val in enumerate(value):
            value_list = [key,val]
            output_cps = os.path.splitext(input_cps)[0] + str(i) + '.cps'
            output_cps = os.path.join(dirname, output_cps)
            copasi_outfile_num = os.path.splitext(copasi_outfile)[0] + str(i) + '.csv'
            change_list += Change_CPS(root, value_list, copasi_outfile=copasi_outfile_num)
            ET.ElementTree(root).write(output_cps, encoding='utf-8')
            root = ET.fromstring(xmlstring)


    return change_list


def Sensitivity_Setup(input_cps, Species_List = ['Nrg'],  ModelParameterGroup = 'Species', jsonfile=''):
    """
    This function sets up global and local sensitivity
    runs using a specified range of perturbations
    """
    #Species_List = ['Nrg','pSTATc','SOCS','pIFNRJ2_pSTATc']
    change_list = []
    copasi_outfile='Time_Course.csv'

    with open(input_cps, 'r') as fp: xmlstring = fp.read()
    xmlstring = re.sub(' xmlns="[^"]+"','',xmlstring,count=1)
    root = ET.fromstring(xmlstring)

    change_list = []
    copasi_outfile='Time_Course.csv'

    output_dirname='SensDir'

    All_Species = GetCpsvalues(input_cps, ModelParameterGroup = ModelParameterGroup)

    Non_Zero_Species = [elem[0] for elem in All_Species if float(elem[1]) > 0.0]

    if Species_List != 'All':
        #Determine if the name of all species matches in Copasi file...
        unmatched_set = [scan_item for scan_item in Species_List if scan_item not in 
                          [species[0] for species in All_Species]]

        #Raise exception otherwise
        if unmatched_set != []:
            raise Exception('Following species/parameters are not found in COPASI file: \n' + '\n'.join(unmatched_set))
        
        #Skip all species in the provided list which have zero initial concentration
        zero_set = [species for species in Species_List if species not in Non_Zero_Species]
        print "Following species concentrations/parameters have zero initial values, hence ignored:\n",'\n'.join(zero_set)
        nonzero_set = list(set(Species_List) - set(zero_set))
    else:
        nonzero_set = list(Non_Zero_Species)

    param_dict = {}

    dirname = output_dirname + ModelParameterGroup
    status = makedirectory(dirname=dirname)
    #Get values and corresponding compartments in the Copasi file
    for index in [ii for ii,i in enumerate(All_Species) if i[0] in nonzero_set and float(i[1]) > 0.0]:
        output_cps = os.path.splitext(input_cps)[0] + All_Species[index][0] + 'Lower' + '.cps'
        output_cps = os.path.join(dirname, output_cps)
        copasi_outfile_num = os.path.splitext(copasi_outfile)[0] + All_Species[index][0] + 'Lower' + '.csv'
        value_list = [All_Species[index][0], 0.9*float(All_Species[index][1])]
        compartment_list = All_Species[index][::2] if ModelParameterGroup == 'Species' else []
        change_list += Change_CPS(root, value_list, compartments=compartment_list, copasi_outfile=copasi_outfile_num)
        ET.ElementTree(root).write(output_cps, encoding='utf-8')
        root = ET.fromstring(xmlstring)
        #change_list += Change_CPS(input_cps, value_list,compartments=compartment_list,output_dirname=dirname, output_cps=output_cps,copasi_outfile=copasi_outfile_num)

        output_cps = os.path.splitext(input_cps)[0] + All_Species[index][0] + 'Higher' + '.cps'
        output_cps = os.path.join(dirname, output_cps)
        copasi_outfile_num = os.path.splitext(copasi_outfile)[0] + All_Species[index][0] + 'Higher' + '.csv'
        value_list = [All_Species[index][0], 1.1*float(All_Species[index][1])]
        compartment_list = All_Species[index][::2] if ModelParameterGroup == 'Species' else []
        change_list += Change_CPS(root, value_list,compartments=compartment_list, copasi_outfile=copasi_outfile_num)
        ET.ElementTree(root).write(output_cps, encoding='utf-8')
        root = ET.fromstring(xmlstring)
        #change_list += Change_CPS(input_cps, value_list,compartments=compartment_list,output_dirname=dirname, output_cps=output_cps,copasi_outfile=copasi_outfile_num)
        param_dict.update({All_Species[index][0]:All_Species[index][1]})

    if jsonfile != '':
        with open(jsonfile,'w') as fp: json.dump(param_dict, fp, indent=4)

    return change_list

def Batch_Run(TargetDir,logfile='_out.log',errfile='_err.log'):
    """
    Checks if directory exists and it not empty
    Then sets up batch run over all COPASI files
    """
    if list(os.walk(TargetDir)) == []: raise Exception('Directory does not exists')
    if list(os.walk(TargetDir))[0][-1] == []: raise Exception('Directory empty')

    filelist = [files for files in list(os.walk(TargetDir))[0][-1] if os.path.splitext(files)[-1] == '.cps']

    for files in sorted(filelist): 
        out, err = OdeRun(os.path.join(TargetDir,files))
        print 'Run completed for file: {0} in directory {1}'.format(files,TargetDir)
        with open(os.path.join(TargetDir, os.path.splitext(files)[0] + logfile),'w') as fp: fp.write(out)
        with open(os.path.join(TargetDir, os.path.splitext(files)[0] + errfile),'w') as fp: fp.write(err)

def Calculate_Sensitivity(TargetDir, input_json):
    """
    Calculates the normalized sensitivity using the data
    obtained from separate simulations
    """
    input_filelist = glob.glob(os.path.join(TargetDir, '*Lower.csv'))
    with open(input_filelist[0], 'r') as fp: header = fp.read().splitlines()[0]
    
    with open(input_json,'r') as fp: Parameter_Value = json.load(fp)

    for ii, files in enumerate(input_filelist):
        m = re.search('\w+Course(\w+)Lower\.csv', files)
        input_parameter = m.group(1)
        input_value = float(Parameter_Value[input_parameter])
        data_slice_lower = np.genfromtxt(files,delimiter=',',skip_header=1)
        files_higher = files.split('Lower.csv')[0] + 'Higher.csv'
        data_slice_upper = np.genfromtxt(files_higher,delimiter=',',skip_header=1)
        
        sens_data = np.abs(data_slice_upper - data_slice_lower)/(0.2*input_value)
        sens_data[:, 0] = data_slice_lower[:, 0]
        format_string = ','.join(['%.1f'] + ['%f' for _ in range(sens_data.shape[1] - 1)])
        np.savetxt(files.split('Lower.csv')[0] + 'Sens.csv', sens_data, fmt=format_string, delimiter=',',header=header)
        print 'Calculation complete for parameter {0}'.format(input_parameter)
#        if ii == 0:
#            data = np.dstack((sens_data,))
#        else:
#            data = np.dstack((data, sens_data))







if __name__ == '__main__':

    copasifile = 'AR_Base_Interfaces_Particle.cps'
    #Setup parameter scans
    change_list = Scan_Setup(copasifile)
    #provide the list of directories
    dirlist = []
    #get all the directories in current folder
    dirlist = [elem[:-1] for elem in glob.glob('ScanDir*/')]

    for dirs in dirlist: Batch_Run(dirs)


