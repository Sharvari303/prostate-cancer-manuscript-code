%% 2020/07
%%%    two cancer cell populations growth model - PRCA draft shared by Ravi
%%% version 2 - update testosterone-dependent growth rate r=r0+pt to hill function

function write_model_v02()

fid             = fopen('merged_model_v02.m','w');

fprintf(fid, '%s\n', 'function dydt = merged_model_v02(t,y,parameters)');

fprintf(fid, '%s\n\n', '');

fprintf(fid, '%s\n', '%% preallocate');
fprintf(fid, '%s\n\n', 'dydt                           = zeros(2,1);');

indices = fileread('define_indices_v02.m');
fprintf(fid, '%s', indices);

parameters = fileread('unpack_parameters_v02.m');
fprintf(fid, '%s', parameters);

rates = fileread('growth_rate_v02.m');
fprintf(fid, '%s', rates);

fclose(fid);
