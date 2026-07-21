/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
#include <algorithm>  // For std::sort

void create_cell_types( void )
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.parameters.o2_proliferation_saturation = 7; //above 7, max rate
	cell_defaults.parameters.o2_proliferation_threshold= 5;   //below 5, zero rate
	cell_defaults.parameters.o2_reference = 38.0; //38.0;  
	cell_defaults.parameters.o2_necrosis_max =2.5 ;
	cell_defaults.parameters.o2_necrosis_threshold = 5;  // 
	cell_defaults.parameters.max_necrosis_rate=0.0028;  //this rate only gets used when stochastic necrosis is active
	cell_defaults.parameters.necrosis_type = PhysiCell_constants::deterministic_necrosis;
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
    Cell_Definition* pCD = find_cell_definition( "PTEN_normal" ); 
	pCD->functions.update_phenotype = PTEN_normal_cellgrowth;

    pCD = find_cell_definition( "PTEN_deleted" ); 
	pCD->functions.update_phenotype = PTEN_deleted_cellgrowth;  
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
		    Cell* pC = NULL;
	    double tumor_radius = parameters.doubles("tumor_radius");
        double cell_radius = cell_defaults.phenotype.geometry.radius; 
		Cell_Definition* pCD = cell_definitions_by_index[0]; //0 stands for S cancer cells
		Cell_Definition* pCD_R = cell_definitions_by_index[1]; //1 stands for R cancer cells
        
		//creating array to store all positions generated 
		std::vector<std::vector<double>> store_positions_S;
		while (store_positions_S.size() < parameters.ints("number_of_cells_S") ) 
        {   //printf("store_positions=%f\n" ,store_positions);
            std::vector<double> position = {0,0,0}; 
			double theta = 2*3.14*UniformRandom();
			position[0] = tumor_radius*UniformRandom()*cos(theta); 
			position[1] = tumor_radius*UniformRandom()*sin(theta);
			position[2] = 0; 
			//printf("position[0]=%f\n" ,position[0]);
		    std::vector<double> compare_position = {0,0,0};  //temporary variable to store position of circle for comparison
			double distance = 0;
			bool overlapping = false;
			for (int j=0 ; j < store_positions_S.size() ; j++) 
			{
               compare_position = store_positions_S[j];
			   distance = distanceCalculate(compare_position[0], position[0], compare_position[1], position[1]);
               
			   if ( distance < 2.0*cell_radius ) 
			   {  //printf("distance=%f\n" ,distance);
			      //printf("2.0*cell_radius=%f\n" ,2.0*cell_radius);
                  overlapping = true;
				  break;
			   }

			}

			if (!overlapping )
			{  //printf("distance=%f\n" ,distance);
			   //printf("2.0*cell_radius=%f\n" ,2.0*cell_radius);
				pC = create_cell( *pCD ); 
			    pC->assign_position( position );
            
			    store_positions_S.push_back(position);
			}

			
		
		}
	    std::vector<std::vector<double>> store_positions_R;
		while (store_positions_R.size() < parameters.ints("number_of_cells_R") ) 
        {   //printf("store_positions=%f\n" ,store_positions);
            std::vector<double> position = {0,0,0}; 
			double theta = 2*3.14*UniformRandom();
			position[0] = tumor_radius*UniformRandom()*cos(theta); 
			position[1] = tumor_radius*UniformRandom()*sin(theta);
			position[2] = 0; 
			//printf("position[0]=%f\n" ,position[0]);
		    std::vector<double> compare_position = {0,0,0};  //temporary variable to store position of circle for comparison
			double distance = 0;
			bool overlapping = false;
			for (int j=0 ; j < store_positions_R.size() ; j++) 
			{
               compare_position = store_positions_R[j];
			   distance = distanceCalculate(compare_position[0], position[0], compare_position[1], position[1]);
               
			   if ( distance < 2.0*cell_radius ) 
			   { // printf("distance=%f\n" ,distance);
			     // printf("2.0*cell_radius=%f\n" ,2.0*cell_radius);
                  overlapping = true;
				  break;
			   }

			}

			if (!overlapping )
			{  //printf("distance=%f\n" ,distance);
			   //printf("2.0*cell_radius=%f\n" ,2.0*cell_radius);
				pC = create_cell( *pCD_R ); 
			    pC->assign_position( position );
            
			    store_positions_R.push_back(position);
			}

			
		
		}
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();
	set_parameters_from_distributions();
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

void PTEN_normal_cellgrowth(Cell* pCell, Phenotype& phenotype, double dt ) //model that determines cell growth rate based on local EGF and testosterone levels
{
    // Get cohort and construct prefix
    int cohort = parameters.ints("cohort");
    std::string prefix;
    if (cohort == 0) prefix = "BR_";
    else if (cohort == 1) prefix = "TR_";
    else if (cohort == 2) prefix = "CTR_";
    else { std::cerr << "Invalid cohort: " << cohort << std::endl; return; }

    // Get testosterone
    static int testosterone_level_index = microenvironment.find_density_index("testosterone");
    double testosterone_level = pCell->nearest_density_vector()[testosterone_level_index];

    // Get parameters with prefix (sensitive parameters)
    double r_s = parameters.doubles(prefix + "rs");
    double map_k = parameters.doubles(prefix + "map_k");
    double a_s = parameters.doubles(prefix + "as");
    double K = parameters.doubles("K");  // Shared
    double ms = parameters.doubles(prefix + "ms");
    double ns = parameters.doubles(prefix + "ns");
    double ps = parameters.doubles(prefix + "ps");


double y_R = 0; //# of PTEN_deleted cells in simulation box
double y_R_density = 0 ; //density of PTEN_deleted cells in simulation box
Cell* pCell_counter = NULL;
static Cell_Definition* pCC = find_cell_definition( "PTEN_deleted" );

// loop over all cells
for ( int i=0; i < (*all_cells).size(); i++ )
{
    pCell_counter = (*all_cells)[i];
    if (pCell_counter->phenotype.death.dead == false && pCell_counter->type == pCC->type)
    {
       y_R++;
	 
    }

}

int mass_one_cell_index = parameters.doubles.find_index("mass_one_cell");
double mass_one_cell = parameters.doubles(mass_one_cell_index);
int tumor_radius_index = parameters.doubles.find_index("tumor_radius");
double tumor_radius = parameters.doubles(tumor_radius_index);

double vol =  3.14*pow(tumor_radius,2)*20*pow(10,-9)*pow(10,6) ; //replace with box dimension inputs 
y_R_density = y_R * mass_one_cell/vol;

//adding effect of oxygen concentration on growth rate 
static int oxygen_substrate_index = pCell->get_microenvironment()->find_density_index( "oxygen" ); 
// sample the microenvironment to get the pO2 value 	
double pO2 = (pCell->nearest_density_vector())[oxygen_substrate_index]; // PhysiCell_constants::oxygen_index]; 

// this multiplier is for linear interpolation of the oxygen value 
	double multiplier = 1.0;
	if( pO2 < pCell->parameters.o2_proliferation_saturation )
	{
		multiplier = ( pO2 - pCell->parameters.o2_proliferation_threshold ) 
			/ ( pCell->parameters.o2_proliferation_saturation - pCell->parameters.o2_proliferation_threshold );
	}
	if( pO2 < pCell->parameters.o2_proliferation_threshold )
	{ 
		multiplier = 0.0; 
	}

pCell->phenotype.cycle.data.transition_rate(0,0)=  multiplier * ( r_s*(1 + (ms/r_s * pow(testosterone_level,ns)) / (pow(ps,ns) + pow(testosterone_level,ns)))*map_k/(30.0*24.0*60.0)); // + (a_s/(30.0*24.0*60.0))*(K - y_R_density)/K );
//printf("sensitive= %f\n" , pCell->phenotype.cycle.data.transition_rate(0,0));
// if dead, set secretion / uptake rates to zero
// trick: if dead, overwrite with NULL function pointer. 

//printf("term_R= %f\n" , (K - y_R_density)/K);

// Update necrosis rate
	multiplier = 0.0;
	if( pO2 < pCell->parameters.o2_necrosis_threshold )
	{
		multiplier = ( pCell->parameters.o2_necrosis_threshold - pO2 ) 
			/ ( pCell->parameters.o2_necrosis_threshold - pCell->parameters.o2_necrosis_max );
	}
	if( pO2 < pCell->parameters.o2_necrosis_max )
	{ 
		multiplier = 1.0; 
	}	
	
	// now, update the necrosis rate 
	int necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model );
	pCell->phenotype.death.rates[necrosis_index] = multiplier * pCell->parameters.max_necrosis_rate; 
	
	// check for deterministic necrosis 
	
	if( pCell->parameters.necrosis_type == PhysiCell_constants::deterministic_necrosis && multiplier > 1e-16 )
	{ pCell->phenotype.death.rates[necrosis_index] = 9e99 ; } 

if( phenotype.death.dead == true )
{
phenotype.secretion.set_all_secretion_to_zero(); 
phenotype.secretion.set_all_uptake_to_zero(); 
pCell->functions.update_phenotype = NULL;
}

return; 
}  

void PTEN_deleted_cellgrowth(Cell* pCell, Phenotype& phenotype, double dt )  //model that determines cell growth rate based on local EGF and testosterone levels
{

	// Get cohort and construct prefix
    int cohort = parameters.ints("cohort");
    std::string prefix;
    if (cohort == 0) prefix = "BR_";
    else if (cohort == 1) prefix = "TR_";
    else if (cohort == 2) prefix = "CTR_";
    else { std::cerr << "Invalid cohort: " << cohort << std::endl; return; }

    // Get testosterone
    static int testosterone_level_index = microenvironment.find_density_index("testosterone");
    double testosterone_level = pCell->nearest_density_vector()[testosterone_level_index];

    // Get parameters with prefix (resistant parameters)
    double r_r = parameters.doubles(prefix + "rr");
    double map_k = parameters.doubles(prefix + "map_k");
    double a_r = parameters.doubles(prefix + "ar");
    double K = parameters.doubles("K");  // Shared
    double mr = parameters.doubles(prefix + "mr");
    double nr = parameters.doubles(prefix + "nr");
    double pr = parameters.doubles(prefix + "pr");

	double y_S = 0; //# of PTEN_normal cells in simulation box
	double y_S_density = 0 ; //density of PTEN_normal cells in simulation box
	Cell* pCell_counter = NULL;
	static Cell_Definition* pCC = find_cell_definition( "PTEN_normal" );

	// loop over all cell
	for ( int i=0; i < (*all_cells).size(); i++ )
	{
		pCell_counter = (*all_cells)[i];
		if (pCell_counter->phenotype.death.dead == false && pCell_counter->type == pCC->type)
		{
		y_S++;
		}

	}

	int mass_one_cell_index = parameters.doubles.find_index("mass_one_cell");
	double mass_one_cell = parameters.doubles(mass_one_cell_index);

	int tumor_radius_index = parameters.doubles.find_index("tumor_radius");
	double tumor_radius = parameters.doubles(tumor_radius_index);

	double vol =  3.14*pow(tumor_radius,2)*20*pow(10,-9)*pow(10,6) ; //replace with box dimension inputs 
	y_S_density = y_S * mass_one_cell/vol;

	//adding effect of oxygen concentration on growth rate 
	static int oxygen_substrate_index = pCell->get_microenvironment()->find_density_index( "oxygen" ); 
	// sample the microenvironment to get the pO2 value 	
	double pO2 = (pCell->nearest_density_vector())[oxygen_substrate_index]; // PhysiCell_constants::oxygen_index]; 

	// this multiplier is for linear interpolation of the oxygen value 
		double multiplier = 1.0;
		if( pO2 < pCell->parameters.o2_proliferation_saturation )
		{
			multiplier = ( pO2 - pCell->parameters.o2_proliferation_threshold ) 
				/ ( pCell->parameters.o2_proliferation_saturation - pCell->parameters.o2_proliferation_threshold );
		}
		if( pO2 < pCell->parameters.o2_proliferation_threshold )
		{ 
			multiplier = 0.0; 
		}

	//printf("pO2= %f\n" , multiplier);

	pCell->phenotype.cycle.data.transition_rate(0,0)= multiplier * (r_r *(1 + (mr/r_r * pow(testosterone_level,nr)) / (pow(pr,nr) + pow(testosterone_level,nr)))*map_k/(30.0*24.0*60.0)) ; //+ (a_r/(30.0*24.0*60.0))*(K - y_S_density)/K );
	//printf("resistant= %f\n" , pCell->phenotype.cycle.data.transition_rate(0,0));
	// if dead, set secretion / uptake rates to zero
	// trick: if dead, overwrite with NULL function pointer 
	//printf("term_S= %f\n" , (K - y_S_density)/K);

	// Update necrosis rate
	multiplier = 0.0;
	if( pO2 < pCell->parameters.o2_necrosis_threshold )
	{
		multiplier = ( pCell->parameters.o2_necrosis_threshold - pO2 )
			/ ( pCell->parameters.o2_necrosis_threshold - pCell->parameters.o2_necrosis_max );
	}
	if( pO2 < pCell->parameters.o2_necrosis_max )
	{
		multiplier = 1.0;
	}

	// now, update the necrosis rate
	int necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model );
	pCell->phenotype.death.rates[necrosis_index] = multiplier * pCell->parameters.max_necrosis_rate;

	// check for deterministic necrosis

	if( pCell->parameters.necrosis_type == PhysiCell_constants::deterministic_necrosis && multiplier > 1e-16 )
	{ pCell->phenotype.death.rates[necrosis_index] = 9e99 ; }

	if(phenotype.death.dead == true)
	{
	phenotype.secretion.set_all_secretion_to_zero(); 
	phenotype.secretion.set_all_uptake_to_zero(); 
	pCell->functions.update_phenotype = NULL;
	}

	return; 
	
}

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)  //different for 2D and 3D cases - this is a 2D case, adapted from an example PhysiCell model 
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*2; //*sqrt(3); (6)
	double y_spacing= cell_radius*2;
	//double z_spacing= cell_radius*sqrt(3);// (1)
	
	std::vector<double> tempPoint(3,0.0);
	// std::vector<double> cylinder_center(3,0.0);
	
	//for(double z=-sphere_radius;z<sphere_radius;z+=z_spacing, zc++) //(2)
	//{ //(3)
		for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
		{
			for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
				tempPoint[1]=y + (xc%2) * cell_radius;
				tempPoint[2]=0; //z; (4)
				
				if(sqrt(norm_squared(tempPoint))< sphere_radius)
				{ cells.push_back(tempPoint); }
			}
			
	    }
	//} //(5)
	return cells; //stores x and y positions of potential centers of circles/agents initialized in the simulation 
	
}

double distanceCalculate(double x1, double y1, double x2, double y2)//calculates distance between centers of two circles/agents 
{
	double x = x1 - x2; //calculating number to square in next step
	double y = y1 - y2;
	double dist;

	dist = pow(x, 2) + pow(y, 2);       //calculating Euclidean distance
	dist = sqrt(dist);                  

	return dist;
}

void combined_analysis(void)
{
    // Get cell definitions
    static Cell_Definition* pTEN_N = find_cell_definition("PTEN_normal");
    static Cell_Definition* pTEN_D = find_cell_definition("PTEN_deleted");

    // Clustering calculations
    double sum_C_all = 0.0, sum_C_N = 0.0, sum_C_D = 0.0;
    int count_all = 0, count_N = 0, count_D = 0;

    // Population counters
    int total_cells = 0;
    int alive_pten_normal = 0;
    int dead_pten_normal = 0;
    int alive_pten_deleted = 0;
    int dead_pten_deleted = 0;

    // Testosterone depletion zone
    double depletion_threshold = 1;  // Adjust if needed
    int testosterone_substrate_index = microenvironment.find_density_index("testosterone");
    double max_radius_testosterone_depleted = 0.0;

    // Front radius: distances from origin (0,0)
    std::vector<double> distances_all, distances_N, distances_D;

    // Loop through all cells
    for (size_t i = 0; i < (*all_cells).size(); ++i)
    {
        Cell* pCell = (*all_cells)[i];
        total_cells++;

        if (pCell->phenotype.death.dead) {
            if (pCell->type == pTEN_N->type) dead_pten_normal++;
            else if (pCell->type == pTEN_D->type) dead_pten_deleted++;
            continue;
        }

        if (pCell->type == pTEN_N->type) alive_pten_normal++;
        else if (pCell->type == pTEN_D->type) alive_pten_deleted++;

        // Calculate radial distance from origin (0,0)
        double radius = sqrt(pCell->position[0] * pCell->position[0] + pCell->position[1] * pCell->position[1]);
        distances_all.push_back(radius);
        if (pCell->type == pTEN_N->type) distances_N.push_back(radius);
        else if (pCell->type == pTEN_D->type) distances_D.push_back(radius);

        // Sample testosterone concentration
        double testosterone_conc = 0.0;
        if (testosterone_substrate_index >= 0) {
            std::vector<double> pos = {pCell->position[0], pCell->position[1], pCell->position[2]};
            testosterone_conc = microenvironment.nearest_density_vector(pos)[testosterone_substrate_index];
        }

        // Update depletion radius
        if (testosterone_conc < depletion_threshold && radius > max_radius_testosterone_depleted) {
            max_radius_testosterone_depleted = radius;
        }

        // Clustering calculations
        int alive_neighbors = 0;
        int same_type_neighbors = 0;
        for (Cell* pNeighbor : pCell->state.neighbors)
        {
            if (pNeighbor->phenotype.death.dead) continue;
            ++alive_neighbors;
            if (pNeighbor->type == pCell->type) ++same_type_neighbors;
        }
        double Ci = 0.0;
        if (alive_neighbors > 0)
        {
            Ci = static_cast<double>(same_type_neighbors) / static_cast<double>(alive_neighbors);
        }
        sum_C_all += Ci;
        ++count_all;
        if (pCell->type == pTEN_N->type) { sum_C_N += Ci; ++count_N; }
        else if (pCell->type == pTEN_D->type) { sum_C_D += Ci; ++count_D; }
    }

    // Compute 95th percentile front radii
    auto percentile_95 = [](std::vector<double>& dists) -> double {
        if (dists.empty()) return 0.0;
        std::sort(dists.begin(), dists.end());
        size_t idx = static_cast<size_t>(0.95 * dists.size());
        return dists[idx];
    };
    double front_radius_all = percentile_95(distances_all);
    double front_radius_N = percentile_95(distances_N);
    double front_radius_D = percentile_95(distances_D);

    // Compute front speeds
    static double prev_front_radius_all = 0.0, prev_front_radius_N = 0.0, prev_front_radius_D = 0.0;
    static double prev_time = 0.0;
    double front_speed_all = 0.0, front_speed_N = 0.0, front_speed_D = 0.0;
    double current_time = PhysiCell_globals.current_time;
    if (prev_time > 0.0) {
        double dt = current_time - prev_time;
        front_speed_all = (front_radius_all - prev_front_radius_all) / dt;
        front_speed_N = (front_radius_N - prev_front_radius_N) / dt;
        front_speed_D = (front_radius_D - prev_front_radius_D) / dt;
    }
    prev_front_radius_all = front_radius_all;
    prev_front_radius_N = front_radius_N;
    prev_front_radius_D = front_radius_D;
    prev_time = current_time;

    // Compute clustering averages
    double C_avg_all = (count_all > 0) ? (sum_C_all / count_all) : 0.0;
    double C_avg_N = (count_N > 0) ? (sum_C_N / count_N) : 0.0;
    double C_avg_D = (count_D > 0) ? (sum_C_D / count_D) : 0.0;

    // Get current time
    double t_min = PhysiCell_globals.current_time;

    // CSV writing
    bool need_header = false;
    {
        std::ifstream test("output/analysis_over_time.csv");
        if (!test.good()) { need_header = true; }
    }
    std::ofstream file("output/analysis_over_time.csv", std::ios::app);
    if (need_header)
    {
        file << "time_min,C_avg_all,C_avg_PTEN_normal,C_avg_PTEN_deleted,";
        file << "Total_Cells,Alive_PTEN_Normal,Dead_PTEN_Normal,Alive_PTEN_Deleted,Dead_PTEN_Deleted,";
        file << "Testosterone_Depletion_Radius,";
        file << "Front_Radius_All,Front_Radius_PTEN_Normal,Front_Radius_PTEN_Deleted,";
        file << "Front_Speed_All,Front_Speed_PTEN_Normal,Front_Speed_PTEN_Deleted\n";
    }
    file << t_min << "," << C_avg_all << "," << C_avg_N << "," << C_avg_D << ",";
    file << total_cells << "," << alive_pten_normal << "," << dead_pten_normal << ",";
    file << alive_pten_deleted << "," << dead_pten_deleted << ",";
    file << max_radius_testosterone_depleted << ",";
    file << front_radius_all << "," << front_radius_N << "," << front_radius_D << ",";
    file << front_speed_all << "," << front_speed_N << "," << front_speed_D << "\n";
    return;
}

