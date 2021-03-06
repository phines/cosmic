%% simulate 39-bus case
clear all; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
end

% simulation time
t_max = 400;

% select data case to simulate
ps = updateps(case39_ps);

% initialize the case
opt = psoptions;
opt.sim.integration_scheme = 1;
opt.sim.dt_default = 0.1;
opt.nr.use_fsolve = true;
opt.verbose = true;
ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);

% build the machine variables
[ps.mac,ps.exc,ps.gov]      = get_mac_state(ps,'salient');
% initialize relays
ps.relay                    = get_relays(ps,'all',opt);

%% build an event matrix
event = zeros(3,C.ev.cols);
% start
event(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
% trip a branch
event(2,[C.ev.time C.ev.type]) = [50 C.ev.trip_branch];
event(2,C.ev.branch_loc) = 32;
% trip a branch
event(3,[C.ev.time C.ev.type]) = [100 C.ev.trip_branch];
event(3,C.ev.branch_loc) = 33;
% trip a branch
event(4,[C.ev.time C.ev.type]) = [200 C.ev.trip_branch];
event(4,C.ev.branch_loc) = 24;
% trip a branch
event(5,[C.ev.time C.ev.type]) = [300 C.ev.trip_branch];
event(5,C.ev.branch_loc) = 23;
% set the end time
event(6,[C.ev.time C.ev.type]) = [t_max C.ev.finish];

%% run the simulation
[outputs,ps] = simgrid(ps,event,'sim_case39',opt);

%% print the results
fname = outputs.outfilename;
[t,delta,omega,Pm,Eap,temp,Vmag,theta] = read_outfile(fname,ps);
omega_0 = 2*pi*ps.frequency;
omega_pu = omega / omega_0;

figure(1); clf; hold on; 
nl = size(omega_pu,2); colorset = varycolor(nl);
set(gca,'ColorOrder',colorset,'FontSize',18,'Xtick',[0 100 200 300 400],...
    'Xlim',[0 410],'Ylim',[0.985 1.07]);
plot(t',omega_pu');
ylabel('\omega (pu)','FontSize',18);
xlabel('time (sec.)','FontSize',18);
print -depsc2 -r600 ~/Desktop/case39_omegapu
%legend(cellstr(num2str((1:nl)', 'omega_%d'))); legend boxoff; 

figure(2); clf; hold on; 
nl = size(theta,2); colorset = varycolor(nl);
set(gca,'ColorOrder',colorset,'FontSize',18,'Xtick',[0 100 200 300 400],...
    'Xlim',[0 410],'Ylim',[-0.30 0.3]);
plot(t',theta');
ylabel('\theta','FontSize',18);
xlabel('time (sec.)','FontSize',18);
print -depsc2 -r600 ~/Desktop/case39_theta
%legend(cellstr(num2str((1:nl)', 'theta_%d'))); legend boxoff;

figure(3); clf; hold on; 
nl = size(Vmag,2); colorset = varycolor(nl);
set(gca,'ColorOrder',colorset,'FontSize',18,'Xtick',[0 100 200 300 400],...
    'Xlim',[0 410],'Ylim',[0.85 1.18]);
plot(t',Vmag');
ylabel('|V|','FontSize',18);
xlabel('time (sec.)','FontSize',18);
print -depsc2 -r600 ~/Desktop/case39_vmag
%legend(cellstr(num2str((1:nl)', 'Vmag_%d'))); legend boxoff;

figure(4); clf; hold on; 
nl = size(temp,2); colorset = varycolor(nl);
set(gca,'ColorOrder',colorset,'FontSize',18,'Xtick',[0 100 200 300 400],...
    'Xlim',[0 410],'Ylim',[0 2.5e4]);
plot(t',temp');
ylabel('temperature','FontSize',18);
xlabel('time (sec.)','FontSize',18);
print -depsc2 -r600 ~/Desktop/case39_temp
%legend(cellstr(num2str((1:nl)', 'temp_%d'))); legend boxoff;
