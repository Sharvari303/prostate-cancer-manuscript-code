%% 2020/07
%%%    two cancer cell populations growth model - PRCA draft shared by Ravi

close all; clear; clc;

%% model simulation

% write_model_v02

load('hiss_fit_BR_EGF10nM.mat')
% pten                = [561000 56100 44880 28050 5610];

% simulation settings
sim_start                   = 0; % month
sim_end                     = 39; % month; %init 0.001: 15; % init 1e-6:65-Sen similar range; 50-Res similar range
sim_increment               = 0.2; % month
tspan                       = sim_start:sim_increment:sim_end;

options                     = odeset('AbsTol',1e-07, 'RelTol',1e-06);


% model set up
parameters_orig             = define_parameters_v02;


p_index_ar                  = 5;
p_index_T                   = 6;
% p_index_T0                  = 4; % for no treatment - when testosterone = 0
p_index_r0r                 = 8;
p_index_mr                  = 9;
p_index_nr                  = 10;
p_index_pr                  = 11;


testconc                    = [1:8];

%% PTEN DEL: 5610 vs normal

parameters                  = parameters_orig;
y_sim                        = NaN(length(tspan),2,length(testconc));

for i = 1:length(testconc)
    parameters(p_index_T)       = testconc(i); % different testosterone concentration

%     if (testconc(i) == 0)
%         parameters(p_index_T)   = parameters(p_index_T0);
%     end
    
    yinit                       = get_init_conditions_v02(parameters);

    [t_sim, y_sim(:,:,i)]              = ode15s(@(t,y) merged_model_v02(t, y, parameters), tspan, yinit, options);

end



%% PTEN DEL: 28050 vs normal

parameters                  = parameters_orig;
parameters(p_index_ar)    	= 0.0084;

parameters(p_index_r0r)    	= R(4,1);
parameters(p_index_mr)     	= R(4,2);
parameters(p_index_nr)    	= R(4,3);
parameters(p_index_pr)    	= R(4,4);

y_sim2                        = NaN(length(tspan),2,length(testconc));

for i = 1:length(testconc)
    parameters(p_index_T)       = testconc(i); % different testosterone concentration

%     if (testconc(i) == 0)
%         parameters(p_index_T)   = parameters(p_index_T0);
%     end
    
    yinit                       = get_init_conditions_v02(parameters);

    [t_sim, y_sim2(:,:,i)]              = ode15s(@(t,y) merged_model_v02(t, y, parameters), tspan, yinit, options);

end


%% PTEN DEL: 44880 vs normal

parameters                  = parameters_orig;
parameters(p_index_ar)     	= 0.02464;

parameters(p_index_r0r)    	= R(3,1);
parameters(p_index_mr)     	= R(3,2);
parameters(p_index_nr)    	= R(3,3);
parameters(p_index_pr)    	= R(3,4);

y_sim3                        = NaN(length(tspan),2,length(testconc));

for i = 1:length(testconc)
    parameters(p_index_T)       = testconc(i); % different testosterone concentration

%     if (testconc(i) == 0)
%         parameters(p_index_T)   = parameters(p_index_T0);
%     end
    
    yinit                       = get_init_conditions_v02(parameters);

    [t_sim, y_sim3(:,:,i)]              = ode15s(@(t,y) merged_model_v02(t, y, parameters), tspan, yinit, options);

end



%% Results

define_VarNames_v01;

Scoeff                                  = 0.145*30; % production rate of R: 0-0.29 ag/mL/cell/day in J. Morken, 2014: DOI: 10.1158/0008-5472.CAN-13-3162)
Rcoeff                                  = 0.1825*30; % production rate of R: 0.015-0.35 ag/mL/cell/day in J. Morken, 2014: DOI: 10.1158/0008-5472.CAN-13-3162)

%%% PTEN DEL 5610 - time course
fig=figure(1);
for i=1:length(testconc)
    subplot(4,2,i)
    plot(t_sim, squeeze(y_sim(:,1,i)),t_sim, squeeze(y_sim(:,2,i)),...
        t_sim, Scoeff*squeeze(y_sim(:,1,i))+Rcoeff*squeeze(y_sim(:,2,i)),'LineWidth',3)
    title(['Testosterone: ' num2str(testconc(i)) 'ng/mL'],'FontSize',14)
    hold on
end
hold off
legend({VarNames{1},VarNames{2},'Total: both S and R cells'},'FontSize',14)
sgtitle("PTEN normal 56100 + PTEN deletion 5610",'FontSize',22)

% Give common xlabel, ylabel and title to figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Months','FontSize',20);
ylabel(han,'PSA (ng/mL)','FontSize',20);
% title(han,'yourTitle');


%%% PTEN DEL 28050 - time course
fig=figure(2);
for i=1:length(testconc)
    subplot(4,2,i)
    plot(t_sim, squeeze(y_sim2(:,1,i)),t_sim, squeeze(y_sim2(:,2,i)),...
        t_sim, Scoeff*squeeze(y_sim2(:,1,i))+Rcoeff*squeeze(y_sim2(:,2,i)),'LineWidth',3)
    title(['Testosterone: ' num2str(testconc(i)) 'ng/mL'],'FontSize',14)
    hold on
end
hold off
legend({VarNames{1},VarNames{2},'Total: both S and R cells'},'FontSize',14)
sgtitle("PTEN normal 56100 + PTEN deletion 28050",'FontSize',22)

% Give common xlabel, ylabel and title to figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Months','FontSize',20);
ylabel(han,'PSA (ng/mL)','FontSize',20);
% title(han,'yourTitle');


%%% PTEN DEL 44880 - time course
fig=figure(3);
for i=1:length(testconc)
    subplot(4,2,i)
    plot(t_sim, squeeze(y_sim3(:,1,i)),t_sim, squeeze(y_sim3(:,2,i)),...
        t_sim, Scoeff*squeeze(y_sim3(:,1,i))+Rcoeff*squeeze(y_sim3(:,2,i)),'LineWidth',3)
    title(['Testosterone: ' num2str(testconc(i)) 'ng/mL'],'FontSize',14)
    hold on
end
hold off
legend({VarNames{1},VarNames{2},'Total: both S and R cells'},'FontSize',14)
sgtitle("PTEN normal 56100 + PTEN deletion 44880",'FontSize',22)

% Give common xlabel, ylabel and title to figure
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Months','FontSize',20);
ylabel(han,'PSA (ng/mL)','FontSize',20);
% title(han,'yourTitle');


%%% comparison bar chart
figure(4)
subplot(3,1,1)
bar(testconc,[y_sim(end,1,1) y_sim2(end,1,1) y_sim3(end,1,1);...
    y_sim(end,1,2) y_sim2(end,1,2) y_sim3(end,1,2);... 
    y_sim(end,1,3) y_sim2(end,1,3) y_sim3(end,1,3);...
    y_sim(end,1,4) y_sim2(end,1,4) y_sim3(end,1,4);...
    y_sim(end,1,5) y_sim2(end,1,5) y_sim3(end,1,5);...
    y_sim(end,1,6) y_sim2(end,1,6) y_sim3(end,1,6);...
    y_sim(end,1,7) y_sim2(end,1,7) y_sim3(end,1,7);...
    y_sim(end,1,8) y_sim2(end,1,8) y_sim3(end,1,8)],'FaceColor',[0.47,0.67,0.19]) 
title (VarNames{1},'FontSize',14)
ylabel('PSA (ng/mL)','FontSize',14)
legend({'PTEN Normal: 56100'},'FontSize',14)
% xlabel('Testosterone (ng/mL)','FontSize',12)

subplot(3,1,2)
bar(testconc,[y_sim(end,2,1) y_sim2(end,2,1) y_sim3(end,2,1);...
    y_sim(end,2,2) y_sim2(end,2,2) y_sim3(end,2,2);...
    y_sim(end,2,3) y_sim2(end,2,3) y_sim3(end,2,3);...
    y_sim(end,2,4) y_sim2(end,2,4) y_sim3(end,2,4);... 
    y_sim(end,2,5) y_sim2(end,2,5) y_sim3(end,2,5);...
    y_sim(end,2,6) y_sim2(end,2,6) y_sim3(end,2,6);...
    y_sim(end,2,7) y_sim2(end,2,7) y_sim3(end,2,7);...
    y_sim(end,2,8) y_sim2(end,2,8) y_sim3(end,2,8)])
title (VarNames{2},'FontSize',14)
ylabel('PSA (ng/mL)','FontSize',14)
legend({'PTEN DEL: 5610 (0.1x)', 'PTEN DEL: 28050 (0.5x)', 'PTEN DEL: 44880 (0.8x)'},'FontSize',14)
% xlabel('Testosterone (ng/mL)','FontSize',12)

subplot(3,1,3)
bar(testconc,[Scoeff*y_sim(end,1,1)+Rcoeff*y_sim(end,2,1) Scoeff*y_sim2(end,1,1)+Rcoeff*y_sim2(end,2,1) Scoeff*y_sim3(end,1,1)+Rcoeff*y_sim3(end,2,1);...
    Scoeff*y_sim(end,1,2)+Rcoeff*y_sim(end,2,2) Scoeff*y_sim2(end,1,2)+Rcoeff*y_sim2(end,2,2) Scoeff*y_sim3(end,1,2)+Rcoeff*y_sim3(end,2,2);...
    Scoeff*y_sim(end,1,3)+Rcoeff*y_sim(end,2,3) Scoeff*y_sim2(end,1,3)+Rcoeff*y_sim2(end,2,3) Scoeff*y_sim3(end,1,3)+Rcoeff*y_sim3(end,2,3);...
    Scoeff*y_sim(end,1,4)+Rcoeff*y_sim(end,2,4) Scoeff*y_sim2(end,1,4)+Rcoeff*y_sim2(end,2,4) Scoeff*y_sim3(end,1,4)+Rcoeff*y_sim3(end,2,4);...
    Scoeff*y_sim(end,1,5)+Rcoeff*y_sim(end,2,5) Scoeff*y_sim2(end,1,5)+Rcoeff*y_sim2(end,2,5) Scoeff*y_sim3(end,1,5)+Rcoeff*y_sim3(end,2,5);...
    Scoeff*y_sim(end,1,6)+Rcoeff*y_sim(end,2,6) Scoeff*y_sim2(end,1,6)+Rcoeff*y_sim2(end,2,6) Scoeff*y_sim3(end,1,6)+Rcoeff*y_sim3(end,2,6);...
    Scoeff*y_sim(end,1,7)+Rcoeff*y_sim(end,2,7) Scoeff*y_sim2(end,1,7)+Rcoeff*y_sim2(end,2,7) Scoeff*y_sim3(end,1,7)+Rcoeff*y_sim3(end,2,7);...
    Scoeff*y_sim(end,1,8)+Rcoeff*y_sim(end,2,8) Scoeff*y_sim2(end,1,8)+Rcoeff*y_sim2(end,2,8) Scoeff*y_sim3(end,1,8)+Rcoeff*y_sim3(end,2,8)])
title ('Total: both S and R cells','FontSize',14)
ylabel('PSA (ng/mL)','FontSize',14)
xlabel('Testosterone (ng/mL)','FontSize',15)
legend({'PTEN DEL: 5610 (0.1x)', 'PTEN DEL: 28050 (0.5x)', 'PTEN DEL: 44880 (0.8x)'},'FontSize',14)

