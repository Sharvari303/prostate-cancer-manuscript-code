%% 2020/07

function yinit = get_init_conditions_v02(parameters)

define_indices_v02;
unpack_parameters_v02

yinit                                   = zeros(2,1);

yinit(index_S)                          =  0.001; % 1e-5; % 0.001;  % ng/ml
yinit(index_R)                          =  0.001; % 1e-5; % 0.001;  % ng/ml


