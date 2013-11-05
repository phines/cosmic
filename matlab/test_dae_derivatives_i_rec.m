%% driver that tests the f and g equations for 9-bus case
clear all; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
end

% select data case to test derivatives
ps = updateps(case9_ps);
%load case2383_mod_ps_dyn;

% set some options
opt = psoptions;
opt.sim.dt_default = 0.1;
opt.nr.use_fsolve = true;
opt.verbose = true;
% initialize the case
ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);

% build the machine variables
[ps.mac,ps.exc,ps.gov] 		= get_mac_state(ps,'salient');
% initialize relays
opt.sim.uvls_limit = 0.9;
ps.relay                    = get_relays(ps,'all',opt);

% build indices
n  = size(ps.bus,1);
ng = size(ps.gen,1);
m  = size(ps.branch,1);
n_sh = size(ps.shunt,1);
ix   = get_indices(n,ng,m,n_sh,opt);

% build xy
[x0,y0] = get_xy_rec(ps,opt);
xy0     = [x0;y0];
xyp0    = zeros(size([x0;y0]));

%% check derivatives for DAE equations
%% check dF_dxy

disp('Checking dF_dxy for x0');
F_handle = @(xy) differential_algebraic_i_rec([],xy,xyp0,ps,opt);
checkDerivatives(F_handle,[],xy0);

xy1 = xy0+randn(size(xy0))*0.1;
disp('Checking dF_dxy for perturbed x0');
F_handle = @(xy) differential_algebraic_i_rec([],xy,xyp0,ps,opt);
checkDerivatives(F_handle,[],xy1);

return
