function [outputs,ps] = sim_case39_n_2(a,b) 
%% simulate 39-bus case
C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end

% simulation time
t_max = 20;

% select data case to simulate
ps = updateps(case39_ps);

% set some options
opt = psoptions;
opt.sim.integration_scheme = 2;
opt.sim.dt_default = 0.1;
opt.nr.use_fsolve = true;
opt.verbose = true;
opt.sim.writelog = true;
% initialize the case
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
event(2,[C.ev.time C.ev.type]) = [10 C.ev.trip_branch];
event(2,C.ev.branch_loc) = a;
% trip a branch
event(3,[C.ev.time C.ev.type]) = [10 C.ev.trip_branch];
event(3,C.ev.branch_loc) = b;
% set the end time
event(4,[C.ev.time C.ev.type]) = [t_max C.ev.finish];

%% run the simulation
[outputs,ps] = simgrid(ps,event,'sim_case39',opt);

%% print the results
%
fname = outputs.outfilename;
[t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd] = read_outfile(fname,ps);
omega_0 = 2*pi*ps.frequency;
omega_pu = omega / omega_0;

figure(1); clf;
subplot(3,1,1); hold on; 
nl = size(omega_pu,2); colorset = varycolor(nl);
set(gca, 'ColorOrder', colorset);
plot(t',omega_pu');
ylabel('omega pu');
%legend(cellstr(num2str((1:nl)', 'omega_%d'))); legend boxoff; 

subplot(3,1,2); hold on; 
nl = size(theta,2); colorset = varycolor(nl);
set(gca, 'ColorOrder', colorset);
plot(t',theta');
ylabel('theta');
%legend(cellstr(num2str((1:nl)', 'theta_%d'))); legend boxoff;

subplot(3,1,3); hold on; 
nl = size(Vmag,2); colorset = varycolor(nl);
set(gca, 'ColorOrder', colorset);
plot(t',Vmag');
ylabel('Vmag');
xlabel('time');
%legend(cellstr(num2str((1:nl)', 'Vmag_%d'))); legend boxoff;

figure(2); clf; 
subplot(3,1,1); hold on; 
nl = size(delta,2); colorset = varycolor(nl);
set(gca, 'ColorOrder', colorset);
plot(t',delta');
ylabel('delta');

subplot(3,1,2); hold on; 
nl = size(Pm,2); colorset = varycolor(nl);
set(gca, 'ColorOrder', colorset);
plot(t',Pm');
ylabel('Pm');

subplot(3,1,3); hold on; 
nl = size(Eap,2); colorset = varycolor(nl);
set(gca, 'ColorOrder', colorset);
plot(t',Eap');
ylabel('Eap');
xlabel('time');

figure(3); clf;
subplot(2,1,1); hold on; 
nl = size(E1,2); colorset = varycolor(nl);
set(gca, 'ColorOrder', colorset);
plot(t',E1');
ylabel('E1');

subplot(2,1,2); hold on; 
nl = size(Efd,2); colorset = varycolor(nl);
set(gca, 'ColorOrder', colorset);
plot(t',Efd');
ylabel('Efd');
%
