%% simulate 9-bus case
clear all; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
end

% select data case to simulate
load case3120sp_ps_dyn;

% initialize the case
opt = psoptions;
opt.nr.use_fsolve = true;
opt.sim.dt_out = 1/30;
opt.verbose = true;
ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);

% build the machine variables
ps.mac  = get_mac_state(ps,'salient');
% initialize relays
ps      = get_relay_state(ps);

%% build an event matrix
event = zeros(3,C.ev.cols);
% start
event(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
% trip a branch
event(2,[C.ev.time C.ev.type]) = [50 C.ev.trip_branch];
event(2,C.ev.branch_loc) = 28;
% set the end time
event(3,[C.ev.time C.ev.type]) = [500 C.ev.finish];

%% run the simulation
[outputs,ps] = simgrid(ps,event,'sim_case3120sp',opt);

%% print the results
fname = outputs.outfilename;
[t,delta,omega,Pm,Eap,temp,Vmag,theta] = read_outfile(fname,ps);
omega_0 = 2*pi*ps.frequency;
omega_pu = omega / omega_0;

figure(1); clf;
subplot(3,1,1); hold on;
plot(t',omega_pu');
ylabel('omega_pu');

subplot(3,1,2);
plot(t',theta');
ylabel('theta');

subplot(3,1,3);
plot(t',Vmag');
ylabel('Vmag');

figure(2); clf;
plot(t',temp');
ylabel('temperature');
