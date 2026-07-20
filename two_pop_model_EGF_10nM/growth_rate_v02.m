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

