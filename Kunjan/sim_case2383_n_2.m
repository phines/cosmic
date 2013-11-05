function [outputs,ps] = sim_case2383_n_2(a,b) 
%% simulate 2383-bus polish case
C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
end

% simulation time
t_max = 3600;

% select data case to simulate
load case2383_mod_ps_dyn;
ps.branch(:,C.br.tap)       = 1;
ps.shunt(:,C.sh.frac_S)     = 0;
ps.shunt(:,C.sh.frac_E)     = 1;
ps.shunt(:,C.sh.frac_Z)     = 0;
ps.shunt(:,C.sh.gamma)      = 0.08;
ps = updateps(ps);


% set some options
opt = psoptions;
opt.sim.integration_scheme = 2;
opt.sim.dt_default = 0.1;
opt.nr.use_fsolve = true;
opt.verbose = true;
opt.sim.writelog = false;
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
[outputs,ps] = simgrid(ps,event,'sim_case2383',opt);

%% print the results
%fname = outputs.outfilename;
%[t,delta,omega,Pm,Eap,temp,Vmag,theta] = read_outfile(fname,ps);
%omega_0 = 2*pi*ps.frequency;
%omega_pu = omega / omega_0;

%figure(1); clf;
%subplot(3,1,1); hold on; 
%nl = size(omega_pu,2); colorset = varycolor(nl);
%set(gca, 'ColorOrder', colorset);
%plot(t',omega_pu');
%ylabel('omega pu');
%legend(cellstr(num2str((1:nl)', 'omega_%d'))); legend boxoff; 

%subplot(3,1,2); hold on; 
%nl = size(theta,2); colorset = varycolor(nl);
%set(gca, 'ColorOrder', colorset);
%plot(t',theta');
%ylabel('theta');
%legend(cellstr(num2str((1:nl)', 'theta_%d'))); legend boxoff;

%subplot(3,1,3); hold on; 
%nl = size(Vmag,2); colorset = varycolor(nl);
%set(gca, 'ColorOrder', colorset);
%plot(t',Vmag');
%ylabel('Vmag');
%xlabel('time');
%legend(cellstr(num2str((1:nl)', 'Vmag_%d'))); legend boxoff;

%figure(2); clf; hold on; 
%nl = size(temp,2); colorset = varycolor(nl);
%set(gca, 'ColorOrder', colorset);
%plot(t',temp');
%ylabel('temperature');
%legend(cellstr(num2str((1:nl)', 'temp_%d'))); legend boxoff;