%% parameters - EGF 100nM

function parameters = define_parameters_v02() 

% Pre-allocate
load('hiss_fit_BR_EGF10nM.mat') % from two_pop_/analysis/pca_ncg_testo_hill.m


parameters                                              = NaN(13,1);
% pten                = [561000 56100 44880 28050 5610];

% Sensitive (S) cell: PTEN normal (56100)
r0s                                      = R(2,1);  
ms                                       = R(2,2);
ns                                       = R(2,3);
ps                                       = R(2,4); 


% Resistent (R) cell: PTEN deletion
a_r                                   	= 0.07; % /month - PTEN: 5610

% T0                                      = 0.133; % g/month - H. Vierhapper, 1997: between[0.09, 0.21], mean=0.133 ± 0.040  (https://doi.org/10.1210/jcem.82.5.3958)
T                                       = 8; % ng/mL, testosterone concentration from leuprolide pk/pd: [min, peak] = [0, 8] ng/ml or [0, 20] nM (nano molar)

K                                       = 200; % ng/mL - local/metastatic values, from PRCA contribution draft

% Scoeff                                  = 0.145*30; % production rate of R: 0-0.29 ag/mL/cell/day in J. Morken, 2014: DOI: 10.1158/0008-5472.CAN-13-3162)
% Rcoeff                                  = 0.1825*30; % production rate of R: 0.015-0.35 ag/mL/cell/day in J. Morken, 2014: DOI: 10.1158/0008-5472.CAN-13-3162)

% resistent (R) cell: PTEN DEL (5610)
r0r                                      = R(5,1);  
mr                                       = R(5,2);
nr                                       = R(5,3);
pr                                       = R(5,4); 

% sensitive cell:
a_s                                       = 0.02576; % /month - PTEN Normal: 56100


% map_k                                    = map_k;

%     parameters(1)                       = r;
    parameters(1)                       = r0s;
    parameters(2)                       = ms;
    parameters(3)                       = ns;
    parameters(4)                       = ps;
    
    parameters(5)                       = a_r;

%     parameters(6)                       = T0;
	parameters(6)                       = T;
    
    parameters(7)                       = K;
     
    parameters(8)                       = r0r;
    parameters(9)                  	= mr;
    parameters(10)                      = nr;
    parameters(11)                    	= pr;

    parameters(12)                    	= map_k;

    parameters(13)                    	= a_s;
   
%     parameters(7)                       = Scoeff;
%     parameters(8)                       = Rcoeff;
    
    
    
   