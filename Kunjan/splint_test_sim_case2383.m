%% simulate 2383-bus polish case
clear all; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
end

addpath ('../numerics');

% simulation time
t_max = 3600;

% select data case to simulate
load case2383_mod_ps_dyn;
ps.branch(:,C.br.tap)       = 1;
ps.shunt(:,C.sh.frac_S)     = 0;
ps.shunt(:,C.sh.frac_E)     = 1;
ps.shunt(:,C.sh.frac_Z)     = 0;
ps.shunt(:,C.sh.gamma)      = 0.08;
ps.branch(10,C.br.status) = 0;
ps.branch(11,C.br.status) = 0;
ps.branch(12,C.br.status) = 0;
ps.branch(113,C.br.status) = 0;
ps.shunt(:,C.sh.P) = 3*ps.shunt(:,C.sh.P);
ps.gen(:,C.ge.P) = 3*ps.gen(:,C.ge.P); 
ps = updateps(ps);

% set the outage to simulate
branch_outages = [15, 106];

% initialize the case
opt = psoptions;
opt.sim.dt_default = 0.1;
opt.sim.writelog = false;
opt.nr.linesearch = 'cubic_spline';
opt.nr.use_fsolve = false;
opt.verbose = true;
 
ps = newpf(ps,opt);
