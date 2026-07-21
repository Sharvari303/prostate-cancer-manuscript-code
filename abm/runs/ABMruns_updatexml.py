#!/usr/bin/env python
# coding: utf-8
# importing the module
import xml.etree.ElementTree as obj
import pandas as pd
import sys

# Hardcoded baseline values for parameters - default starting xml file has baseline values for all parameters of uptake, adh, motility.
BASELINE_UPTAKE_RATE = 1e-05
BASELINE_ADHESION = 0.4
BASELINE_BM_ADHESION = 4 
BASELINE_SPEED = 0.4

filename=sys.argv[1]
excel_sheet=sys.argv[2]
run_number=int(sys.argv[3])

df= pd.read_csv(excel_sheet)

##parsing through the config xml file and updating parameters
tree=obj.parse(filename)
root = tree.getroot()

## Update cohort
cohort_map = {'BR': '0', 'TR': '1', 'CTRL': '2'}
cohort_elem = root.find("./user_parameters/cohort")
if cohort_elem is not None:
    cohort_elem.text = cohort_map.get(df.loc[run_number].at['Cohort'], '0')
    
## Change random_seed
root_random_seed=root.findall("./user_parameters/random_seed")
root_random_seed[0].text=str(df.loc[run_number].at["Seed"])
#print(root_random_seed[0].text)

## Change # of S cells
root_S=root.findall("./user_parameters/number_of_cells_S")
root_S[0].text=str(df.loc[run_number].at["PTEN_normal"]) 
#print(root_S[0].text)

## Change # of R cells
root_R=root.findall("./user_parameters/number_of_cells_R")
root_R[0].text=str(df.loc[run_number].at["PTEN_null"])
#print(root_R[0].text)

# Update uptake rates
multiplier = float(df.loc[run_number].at['Uptake_rate_multiplier']) if df.loc[run_number].at['Uptake_rate_multiplier'] != '-' else 1.0
celltype = df.loc[run_number].at['Uptake_celltype']

if celltype == 'both':
    # Update both
    uptake_normal = root.find("./cell_definitions/cell_definition[@name='PTEN_normal']/phenotype/secretion/substrate[@name='testosterone']/uptake_rate")
    if uptake_normal is not None:
        uptake_normal.text = str(BASELINE_UPTAKE_RATE * multiplier)
    uptake_deleted = root.find("./cell_definitions/cell_definition[@name='PTEN_deleted']/phenotype/secretion/substrate[@name='testosterone']/uptake_rate")
    if uptake_deleted is not None:
        uptake_deleted.text = str(BASELINE_UPTAKE_RATE * multiplier)

elif celltype == 'PTEN_null':
    # Update only deleted
    uptake_deleted = root.find("./cell_definitions/cell_definition[@name='PTEN_deleted']/phenotype/secretion/substrate[@name='testosterone']/uptake_rate")
    if uptake_deleted is not None:
        uptake_deleted.text = str(BASELINE_UPTAKE_RATE * multiplier)

#Update adhesion and motility parameters
adhesion_multiplier = float(df.loc[run_number].at['Cell_cell_adhesion_multiplier']) if df.loc[run_number].at['Cell_cell_adhesion_multiplier'] != '-' else 1.0
speed_multiplier = float(df.loc[run_number].at['Migration_speed_multiplier']) if df.loc[run_number].at['Migration_speed_multiplier'] != '-' else 1.0
adh_celltype = df.loc[run_number].at['AdhMot_celltype']

if adh_celltype == 'PTEN_null':
    
    speed_normal = root.find("./cell_definitions/cell_definition[@name='PTEN_deleted']/phenotype/motility/speed")
    if speed_normal is not None:
        speed_normal.text = str(BASELINE_SPEED * speed_multiplier)

    adhesion_normal = root.find("./cell_definitions/cell_definition[@name='PTEN_deleted']/phenotype/mechanics/cell_cell_adhesion_strength")
    if adhesion_normal is not None:
        adhesion_normal.text = str(BASELINE_ADHESION * adhesion_multiplier)
    
    bm_adhesion_normal = root.find("./cell_definitions/cell_definition[@name='PTEN_deleted']/phenotype/mechanics/cell_BM_adhesion_strength")
    if bm_adhesion_normal is not None:
        bm_adhesion_normal.text = str(BASELINE_BM_ADHESION * adhesion_multiplier)
    

# Update testosterone BC based on Androgen_condition
test_bc_value = '8' if df.loc[run_number].at['Androgen_condition'] == 'High' else '1'
test_var = root.find("./microenvironment_setup/variable[@name='testosterone']")
if test_var is not None:
    # Update initial_condition
    ic = test_var.find("./initial_condition")
    if ic is not None:
        ic.text = test_bc_value
    # Update Dirichlet_boundary_condition
    dbc = test_var.find("./Dirichlet_boundary_condition")
    if dbc is not None:
        dbc.text = test_bc_value
    # Update all boundary_value in Dirichlet_options
    dirichlet_options = test_var.find("./Dirichlet_options")
    if dirichlet_options is not None:
        for bv in dirichlet_options.findall("./boundary_value"):
            bv.text = test_bc_value

tree.write(filename)
