%% simulate 39-bus case
clear all; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end

% simulation time
t_max = 100;

% select data case to simulate
ps = updateps(case39_ps);
% ps = replicate_case(ps,2);          
% ps = unify_generators(ps); 

% set some options
opt = psoptions;
opt.sim.integration_scheme = 1;
opt.sim.dt_default = 1/10;
opt.nr.use_fsolve = true;
% opt.pf.linesearch = 'cubic_spline';
opt.verbose = true;
opt.sim.gen_control = 1;        % 0 = generator without exciter and governor, 1 = generator with exciter and governor
opt.sim.angle_ref = 0;          % 0 = delta_sys, 1 = center of inertia---delta_coi
                                % Center of inertia doesn't work when having islanding
opt.sim.COI_weight = 0;         % 1 = machine inertia, 0 = machine MVA base(Powerworld)
opt.sim.time_delay_ini = 0.5;     % 1 sec delay for each relay. We might set differernt intitial contidtion for different relays in the future.
% Don't forget to change this value (opt.sim.time_delay_ini) in solve_dae.m

% initialize the case
ps = newpf_rec(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
ps = update_load_freq_source(ps);
% build the machine variables
[ps.mac,ps.exc,ps.gov] 		= get_mac_state(ps,'salient');
% initialize relays
ps.relay                    = get_relays(ps,'all',opt);

global t_delay t_prev_check num_ls dist2threshold state_a
n    = size(ps.bus,1);
ng   = size(ps.mac,1);
m    = size(ps.branch,1);
n_sh = size(ps.shunt,1);
ix   = get_indices(n,ng,m,n_sh,opt);
t_delay = inf(size(ps.relay,1),1);
t_delay([ix.re.uvls,ix.re.ufls,ix.re.dist])= opt.sim.time_delay_ini;
t_prev_check = nan(size(ps.relay,1),1);
num_ls = 0;
dist2threshold = inf(size(ix.re.oc,2)*2,1);
state_a = zeros(size(ix.re.oc,2)*2,1);

%% build an event matrix
event = zeros(6,C.ev.cols);
% start
event(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
% trip a branch
event(2,[C.ev.time C.ev.type]) = [1 C.ev.trip_branch];
event(2,C.ev.branch_loc) = 22;
% trip a branch
event(3,[C.ev.time C.ev.type]) = [3 C.ev.trip_branch];
event(3,C.ev.branch_loc) = 33;
% trip a branch
% event(4,[C.ev.time C.ev.type]) = [10 C.ev.trip_branch];
% event(4,C.ev.branch_loc) = 24;
% trip a branch
% event(5,[C.ev.time C.ev.type]) = [3 C.ev.trip_branch];
% event(5,C.ev.branch_loc) = 23;
% % close a branch
% event(3,[C.ev.time C.ev.type]) = [3.1 C.ev.close_branch];
% event(3,C.ev.branch_loc) = 32;
% trip a shunt
% event(2,[C.ev.time C.ev.type]) = [5 C.ev.trip_shunt];
% event(2,C.ev.shunt_loc) = 1;
% close a shunt
% event(3,[C.ev.time C.ev.type]) = [3.1 C.ev.close_shunt];
% event(3,C.ev.shunt_loc) = 1;
% set the end time
event(6,[C.ev.time C.ev.type]) = [t_max C.ev.finish];

%% run the simulation
[outputs,ps] = simgrid_rec(ps,event,'sim_case39',opt);

%% print the results
fname = outputs.outfilename;
[t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd] = read_outfile_rec(fname,ps,opt);
omega_0 = 2*pi*ps.frequency;
omega_pu = omega / omega_0;

CaseName = 'Case39';
Contingency = 'Branch32Tripping';
% Contingency = 'SS';
% Contingency = 'Branch32TrippingReclose';
% Contingency = 'Shunt1TrippingReclose';
% Contingency = 'Shunt1Tripping';
% Contingency = 'LoadSheddingBus5_5P';

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


% figure(4); clf; hold on; 
% nl = size(temp,2); colorset = varycolor(nl);
% % set(gca,'ColorOrder',colorset,'FontSize',18,'Xtick',[0 600 1200 1800],...
% %     'Xlim',[0 50],'Ylim',[0 2000]);
% plot(t,temp);
% ylabel('temperature','FontSize',18);
% xlabel('time (sec.)','FontSize',18);
% print -depsc2 -r600 ~/Desktop/Meetings/case9_temp_R_Branch7_wo_control
% %legend(cellstr(num2str((1:nl)', 'temp_%d'))); legend boxon;

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

% RecResults.Time = t;
% RecResults.Omega_pu = omega_pu;
% RecResults.Theta = theta;
% RecResults.Vmag = Vmag;
% RecResults.Pm = Pm;
% RecResults.Delta = delta;
% RecResults.Eap = Eap;
% RecResults.E1 = E1;
% RecResults.Efd = Efd;
% ResultsStr = sprintf('RecResults_%s_%s_%s_%s_%s_%s.mat',CaseName, Contingency, Control, integration,Ref,Weight);
% save(ResultsStr, '-struct', 'RecResults')
