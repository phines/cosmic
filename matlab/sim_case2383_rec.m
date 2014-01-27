%% simulate 2383-bus polish case
clear all; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end

% simulation time
t_max = 12;

% select data case to simulate
load case2383_mod_ps_dyn;
ps.branch(:,C.br.tap)       = 1;
ps.shunt(:,C.sh.frac_S)     = 0;
ps.shunt(:,C.sh.frac_E)     = 1;
ps.shunt(:,C.sh.frac_Z)     = 0;
ps.shunt(:,C.sh.gamma)      = 0.08;
ps = updateps(ps);

% set the outage to simulate
branch_outages = 15; % [15, 106];

% initialize the case
opt = psoptions2383;
opt.sim.dt_default = 0.1;
opt.sim.writelog = true;
opt.nr.use_fsolve = true;
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
[ps.mac,ps.exc,ps.gov]      = get_mac_state(ps,'salient');
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

%% build an event matrix with branch outages
n_event = length(branch_outages)+2;
event = zeros(n_event,C.ev.cols);
% start
event(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
% trip branches
for i = 1:length(branch_outages)
    event(i+1,[C.ev.time C.ev.type]) = [10 C.ev.trip_branch];
    event(i+1,C.ev.branch_loc) = branch_outages(i);
end
% set the end time
event(end,[C.ev.time C.ev.type]) = [t_max C.ev.finish];

%% run the simulation
[outputs,ps] = simgrid_rec(ps,event,'sim_case2383',opt);

%% print the results
if opt.sim.writelog 
    fname = outputs.outfilename;
    [t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd] = read_outfile_rec(fname,ps);
    omega_0 = 2*pi*ps.frequency;
    omega_pu = omega / omega_0;
    
    figure(1); clf;
    %subplot(3,1,1); hold on;
    nl = size(omega_pu,2); colorset = varycolor(nl);
    set(gca, 'ColorOrder', colorset);
    plot(t',omega_pu');
    ylabel('omega pu');
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
    
end