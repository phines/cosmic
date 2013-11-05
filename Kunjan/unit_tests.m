
%% run unit test that checks results of cosmic R0.9
clear all; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
end

% simulation time
t_max = 1800;

fprintf('----------------------------------------------------\n');
fprintf('COSMIC R0.09 Unit Test \n');

%% unit test for case9
casename = 'case9';
fprintf('----------------------------------------------------\n');
fprintf('Running unit test for %s ... \n',casename);
fprintf('----------------------------------------------------\n');

% select data case to simulate
ps  = updateps(case9_ps);

% read the event data
event = unit_test_case9_events;

% initialize the case
opt = psoptions;
opt.sim.dt_default = 0.1;
opt.nr.use_fsolve = true;
opt.verbose = true;
ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);

% build the machine variables
[ps.mac,ps.exc,ps.gov] = get_mac_state(ps,'salient');
% initialize relays
ps.relay = get_relays(ps);

% run the simulation
[outputs,ps] = simgrid(ps,event,'sim_case9',opt);

% compare the results to those in the file
fname1 = outputs.outfilename;
fname2 = 'sim_case9_unit_test_09.csv';
unit_test_compare(ps,fname1,fname2,casename);

%% do the unit test for case39
casename = 'case39';
fprintf('----------------------------------------------------\n');
fprintf('Running unit test for case39 ... \n');
fprintf('----------------------------------------------------\n');

% select data case to simulate
ps = updateps(case39_ps);

% initialize the case
opt = psoptions;
opt.sim.dt_default = 0.1;
opt.nr.use_fsolve = true;
opt.verbose = true;
ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);

% build the machine variables
[ps.mac,ps.exc,ps.gov] 		= get_mac_state(ps,'salient');
% initialize relays
ps.relay    				= get_relays(ps);

% build an event matrix
event = unit_test_case39_events;

% run the simulation
[outputs,ps] = simgrid(ps,event,'sim_case39',opt);

% compare the results to those in the file
fname1 = outputs.outfilename;
fname2 = 'sim_case39_unit_test_09.csv';
unit_test_compare(ps,fname1,fname2,casename);


