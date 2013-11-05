%% driver that tests the f and g equations for 9-bus case
clear all; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end

% select data case
ps  = updateps(case39_ps);

% initialize the case
opt = psoptions;
opt.nr.use_fsolve = true;
opt.verbose = true;
opt.sim.gen_control = 1;        % 0 = generator without exciter and governor, 1 = generator with exciter and governor
opt.sim.angle_ref = 1;          % 0 = delta_sys, 1 = center of inertia---delta_coi
opt.sim.COI_weight = 0;  


ps = newpf_rec(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);

% build the machine variables
[ps.mac,ps.exc,ps.gov] 		= get_mac_state(ps,'salient');
% initialize relays
ps.relay    				= get_relays(ps,'temperature');

% build indices
n = size(ps.bus,1);
ng = size(ps.gen,1);
m = size(ps.branch,1);
n_sh = size(ps.shunt,1);

% build x and y
[x0,y0] = get_xy_rec(ps,opt);
[x1,y1] = get_xy(ps,opt);

% check derivatives for differential equations
%% check df_dx
x = x0; y = y0;
disp('Checking df_dx_rec for x0');
[f0,df_dx_rec,df_dy_rec] = differential_eqs_rec(0,x,y,ps,opt);
f_handle_rec = @(x) differential_eqs_rec(0,x,y,ps,opt);
checkDerivatives(f_handle_rec,df_dx_rec,x);
% repeat from a perturbed starting point
disp('Checking df_dx_rec for perturbed x0');
x = x0+randn(size(x0))*0.1;
[f0,df_dx_rec,df_dy_rec] = differential_eqs_rec(0,x,y,ps,opt);
checkDerivatives(f_handle_rec,df_dx_rec,x);

% x = x1; y = y1;
% disp('Checking df_dx_polar for x1');
% [f0,df_dx_polar,df_dy_polar] = differential_eqs(0,x,y,ps,opt);
% f_handle_polar = @(x) differential_eqs(0,x,y,ps,opt);
% checkDerivatives(f_handle_polar,df_dx_polar,x);
% % repeat from a perturbed starting point
% disp('Checking df_dx_polar for perturbed x1');
% x = x0+randn(size(x1))*0.1;
% [f0,df_dx_polar,df_dy_polar] = differential_eqs(0,x,y,ps,opt);
% checkDerivatives(f_handle_polar,df_dx_polar,x);


%% check df_dy
x = x0; y = y0;
disp('Checking df_dy_rec for y0');
[f0,df_dx_rec,df_dy_rec] = differential_eqs_rec(0,x,y,ps,opt);
f_handle_rec = @(y) differential_eqs_rec(0,x,y,ps,opt);
checkDerivatives(f_handle_rec,df_dy_rec,y);
% repeat from a perturbed starting point
disp('Checking df_dy_rec for perturbed y0');
y = y0+randn(size(y0))*0.1;
[f0,df_dx_rec,df_dy_rec] = differential_eqs_rec(0,x,y,ps,opt);
checkDerivatives(f_handle_rec,df_dy_rec,y);

% x = x1; y = y1;
% disp('Checking df_dy_polar for y1');
% [f0,df_dx_polar,df_dy_polar] = differential_eqs(0,x,y,ps,opt);
% f_handle_polar = @(y) differential_eqs(0,x,y,ps,opt);
% checkDerivatives(f_handle_polar,df_dy_polar,y);
% % repeat from a perturbed starting point
% disp('Checking df_dy_polar for perturbed y1');
% y = y0+randn(size(y0))*0.1;
% [f0,df_dx_polar,df_dy_polar] = differential_eqs(0,x,y,ps,opt);
% checkDerivatives(f_handle_polar,df_dy_polar,y);


%% check dg_dx
x = x0; y = y0;
disp('Checking dg_dx_rec for x0');
[f0,dg_dx_rec,dg_dy_rec] = algebraic_eqs_rec(0,x,y,ps,opt);
g_handle_rec = @(x) algebraic_eqs_rec(0,x,y,ps,opt);
checkDerivatives(g_handle_rec,dg_dx_rec,x);

% repeat from a perturbed starting point
disp('Checking dg_dx_rec for perturbed x0');
x = x0+randn(size(x0))*0.1;
[f0,dg_dx_rec,dg_dy_rec] = algebraic_eqs_rec(0,x,y,ps,opt);
checkDerivatives(g_handle_rec,dg_dx_rec,x);

% x = x1; y = y1;
% disp('Checking dg_dx_polar for x1');
% [f1,dg_dx_polar,dg_dy_polar] = algebraic_eqs(0,x,y,ps,opt);
% g_handle_polar = @(x) algebraic_eqs(0,x,y,ps,opt);
% checkDerivatives(g_handle_polar,dg_dx_polar,x);
% 
% % repeat from a perturbed starting point
% disp('Checking dg_dx_polar for perturbed x1');
% x = x1+randn(size(x1))*0.1;
% [f1,dg_dx_polar,dg_dy_polar] = algebraic_eqs(0,x,y,ps,opt);
% checkDerivatives(g_handle_polar,dg_dx_polar,x);


%% check dg_dy
x = x0; y = y0;
disp('Checking dg_dy for y0');
[f0,dg_dx_rec,dg_dy_rec] = algebraic_eqs_rec(0,x,y,ps,opt);
g_handle_rec = @(y) algebraic_eqs_rec(0,x,y,ps,opt);
checkDerivatives(g_handle_rec,dg_dy_rec,y);

% repeat from a perturbed starting point
disp('Checking dg_dy for perturbed y0');
y = y0+randn(size(y0))*0.1;
[f0,dg_dx_rec,dg_dy_rec] = algebraic_eqs_rec(0,x,y,ps,opt);
checkDerivatives(g_handle_rec,dg_dy_rec,y);
% 
% x = x1; y = y1;
% disp('Checking dg_dy_polar for y1');
% [f1,dg_dx_polar,dg_dy_polar] = algebraic_eqs(0,x,y,ps,opt);
% g_handle_polar = @(y) algebraic_eqs(0,x,y,ps,opt);
% checkDerivatives(g_handle_polar,dg_dy_polar,y);
% 
% % repeat from a perturbed starting point
% disp('Checking dg_dy_polar for perturbed y1');
% y = y1+randn(size(y1))*0.1;
% [f1,dg_dx_polar,dg_dy_polar] = algebraic_eqs(0,x,y,ps,opt);
% checkDerivatives(g_handle_polar,dg_dy_polar,y);

return
