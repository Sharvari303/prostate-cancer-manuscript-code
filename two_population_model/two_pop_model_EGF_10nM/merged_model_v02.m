function dydt = merged_model_v02(t,y,parameters)


%% preallocate
dydt                           = zeros(2,1);

%% indices

index_S                     = 1;
index_R                     = 2;


%% parameters
r0s                             = parameters(1);
ms                              = parameters(2);
ns                              = parameters(3);
ps                              = parameters(4);
    
a_r                           	= parameters(5);

% T0                            	= parameters(6);
T                              	= parameters(6);

K                              	= parameters(7);
    
r0r                             = parameters(8);
mr                              = parameters(9);
nr                              = parameters(10);
pr                              = parameters(11);


map_k                           = parameters(12);

a_s                             = parameters(13);



%% two cancer cell populations growth model - PRCA draft shared by Ravi
%%% updated with hill function for growth

%%% testosterone-sensitive population S grows with testosterone-dependent growth rate r: 
% In the absence of any therapy, a standard (and constant) T0 value is assumed
% r =             r0 + pT; - probability - v1, updated with hill function in v2


%%% In case of ADT, two scenarios
% LHRH –agonist treatment:
dydt(index_S)           = (map_k *(r0s * (1 + (ms/r0s * T^ns) / (ps^ns + T^ns))) + a_s * (K - y(index_R))/K) * y(index_S);


% LHRH and combined therapy:
dydt(index_R)           = (map_k *(r0r * (1 + (mr/r0r * T^nr) / (pr^nr + T^nr))) + a_r * (K - y(index_S))/K) * y(index_R);

