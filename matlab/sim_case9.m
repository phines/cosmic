%% simulate 9-bus case
clear all; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end

% simulation time
t_max = 30;

% select data case to simulate
ps = updateps(case9_ps);

% set some options
opt = psoptions;
opt.sim.integration_scheme = 1;
opt.sim.dt_default = 1/30;
opt.nr.use_fsolve = true;
% opt.pf.linesearch = 'cubic_spline';
opt.verbose = true;
opt.sim.gen_control = 1;        % 0 = generator without exciter and governor, 1 = generator with exciter and governor
opt.sim.angle_ref = 0;          % 0 = delta_sys, 1 = center of inertia---delta_coi
                                % Now, center of inertia doesn't work when having islanding
opt.sim.COI_weight = 0;         % 1 = machine inertia, 0 = machine MVA base(Powerworld)
opt.sim.time_delay_ini = 0.5;     % 1 sec delay for each relay. We might set differernt intitial contidtion for different relays in the future.

% initialize the case
ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
ps = update_load_freq_source(ps);
% build the machine variables
[ps.mac,ps.exc,ps.gov] 		= get_mac_state(ps,'salient');
% initialize relays
ps.relay                    = get_relays(ps,'all',opt);

global t_delay t_prev_check
t_delay = ones(size(ps.relay,1),1)*opt.sim.time_delay_ini;
t_prev_check = nan(size(ps.relay,1),1);

%% build an event matrix
event = zeros(4,C.ev.cols);
% start
event(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
% trip a branch
% event(2,[C.ev.time C.ev.type]) = [10 C.ev.trip_branch];
% event(2,C.ev.branch_loc) = 2;
% trip a branch
event(2,[C.ev.time C.ev.type]) = [10 C.ev.trip_branch];
event(2,C.ev.branch_loc) = 7;
% % close a branch
% event(3,[C.ev.time C.ev.type]) = [10.1 C.ev.close_branch];
% event(3,C.ev.branch_loc) = 7;
% % Load shedding
% event(2,[C.ev.time C.ev.type]) = [1 C.event.shed_load];
% event(2,C.ev.shunt_loc) = 1;
% event(2,C.event.change_by) = 1;       % 1 = by percentage, 0 = by amount
% event(2,C.ev.quantity) = 1;        % -100%~100% for change by percentage, real MW for change by amount
% Load restoration
% event(3,[C.ev.time C.ev.type]) = [1.1 C.event.shed_load];
% event(3,C.ev.shunt_loc) = 1;
% event(3,C.event.change_by) = 1;       % 1 = by percentage, 0 = by amount
% event(3,C.ev.quantity) = -1;       % -100%~100% for change by percentage, real MW for change by amount

% % trip a shunt
% event(2,[C.ev.time C.ev.type]) = [10 C.ev.trip_shunt];
% event(2,C.ev.shunt_loc) = 1;
% % close a shunt
% event(3,[C.ev.time C.ev.type]) = [10.1 C.ev.close_shunt];
% event(3,C.ev.shunt_loc) = 1;

% set the end time
event(4,[C.ev.time C.ev.type]) = [t_max C.ev.finish];

%% run the simulation
[outputs,ps] = simgrid(ps,event,'sim_case9',opt);


%% print the results
fname = outputs.outfilename;
[t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd] = read_outfile(fname,ps,opt);
omega_0 = 2*pi*ps.frequency;
omega_pu = omega / omega_0;

CaseName = 'Case9';
% Contingency = 'Branch7Tripping';
% Contingency = 'Branch2Tripping';`
% Contingency = 'SS';
% Contingency = 'Branch7TrippingReclose';
% Contingency = 'Branch6TrippingReclose';
% Contingency = 'LoadSheddingBus5_5P';
% Contingency = 'LoadSheddingBus5_100P';
% Contingency = 'LoadShed_Restore_Bus5_5P';
% Contingency = 'LoadSheddingBus5_10P';
% Contingency = 'Branch7-6Trip';


if opt.sim.gen_control; Control = 'w_control'; else Control = 'wo_control'; end
if opt.sim.angle_ref;
    Ref = 'dCOI';
    if opt.sim.COI_weight; Weight = 'Ms'; else Weight = 'MVA'; end
else
    Ref = 'dSYS';
    Weight = 'noWeight';
end


figure(1); clf; hold on; 
nl = size(omega_pu,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xtick',[0 600 1200 1800],...
%     'Xlim',[0 50],'Ylim',[0.995 1.008]);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.995 1.008]);
plot(t,omega_pu);
ylabel('\omega (pu)','FontSize',18);
xlabel('time (sec.)','FontSize',18);
% PrintStr = sprintf('OmegaPu_P_%s_%s_%s',CaseName, Contingency, Control);
% print('-dpng','-r600',PrintStr)

% print -depsc2 -r600 ~/Desktop/Meetings/case9_omegapu_relay_R_Branch7_wo_control
% legend(cellstr(num2str((1:nl)', 'omega_%d'))); legend boxon; 

figure(2); clf; hold on; 
nl = size(theta,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[-0.2 0.5]);
plot(t,theta);
ylabel('\theta','FontSize',18);
xlabel('time (sec.)','FontSize',18);
% PrintStr = sprintf('Theta_P_%s_%s_%s',CaseName, Contingency, Control);
% print('-dpng','-r600',PrintStr)

% print -depsc2 -r600 ~/Desktop/Meetings/case9_theta_R_Branch7_wo_control
% legend(cellstr(num2str((1:nl)', 'theta_%d'))); legend boxon;

figure(3); clf; hold on; 
nl = size(Vmag,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
plot(t,Vmag);
ylabel('|V|','FontSize',18);
xlabel('time (sec.)','FontSize',18);
% PrintStr = sprintf('Vmag_P_%s_%s_%s',CaseName, Contingency, Control);
% print('-dpng','-r600',PrintStr)

% print -depsc2 -r600 ~/Desktop/Meetings/case9_vmag_R_Branch7_wo_control
% legend(cellstr(num2str((1:nl)', 'Vmag_%d'))); legend boxon;

% figure(5); clf; hold on; 
% nl = size(Pm,2); colorset = varycolor(nl);
% % set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
% plot(t,Pm);
% ylabel('Pm','FontSize',18);
% xlabel('time (sec.)','FontSize',18);

figure(6); clf; hold on; 
nl = size(delta,2); colorset = varycolor(nl);
% set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
% plot(t',delta'.*180./pi);
plot(t,delta);
ylabel('Delta','FontSize',18);
xlabel('time (sec.)','FontSize',18);

% figure(7); clf; hold on; 
% nl = size(Eap,2); colorset = varycolor(nl);
% % set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
% plot(t,Eap);
% ylabel('Eap','FontSize',18);
% xlabel('time (sec.)','FontSize',18);
% 
% figure(8); clf; hold on; 
% nl = size(E1,2); colorset = varycolor(nl);
% % set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
% plot(t,E1);
% ylabel('E1','FontSize',18);
% xlabel('time (sec.)','FontSize',18);
% 
% figure(9); clf; hold on; 
% nl = size(Efd,2); colorset = varycolor(nl);
% % set(gca,'ColorOrder',colorset,'FontSize',18,'Xlim',[0 50],'Ylim',[0.88 1.08]);
% plot(t,Efd);
% ylabel('Efd','FontSize',18);
% xlabel('time (sec.)','FontSize',18);

if opt.sim.integration_scheme == 1
    integration = 'Trap';
elseif opt.sim.integration_scheme == 2
    integration = 'ode15i';
elseif opt.sim.integration_scheme == 3
    integration = 'ode15s';
end

% PolarResults.Time = t;
% PolarResults.Omega_pu = omega_pu;
% PolarResults.Theta = theta;
% PolarResults.Vmag = Vmag;
% PolarResults.Temp = temp;
% PolarResults.Pm = Pm;
% PolarResults.Delta = delta;
% PolarResults.Eap = Eap;
% PolarResults.E1 = E1;
% PolarResults.Efd = Efd;
% ResultsStr = sprintf('PolarResults_%s_%s_%s_%s_%s_%s.mat',CaseName, Contingency, Control, integration,Ref,Weight);
% save(ResultsStr, '-struct', 'PolarResults')

