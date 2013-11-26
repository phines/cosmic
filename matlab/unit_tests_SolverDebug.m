% This code is used to debug integration stalling issue within trapezoidal
% rule method.

clear all; close all; clc; C = psconstants;
% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end

% simulation time
t_max = 180;
ps = updateps(case9_ps);
opt = psoptions;

% options under test
opt.sim.integration_scheme = 1;      % 1-Trapezoidal rule; 2-ODE15i; 3-ODE15s 
opt.pf.linesearch = 'cubic_spline';  % default linesearch is 'backtrack'. Options: 'backtrack', 'exact', and 'cubic_spline'
LoadType = 1;                        % 1: exponential load [V^(0.08)]; 0: constant PQ load 
opt.sim.gen_control = 1;             % 0 = generator without exciter and governor, 1 = generator with exciter and governor
check_derivatives = 1;               % 1 = check dae derivatives

% other options
opt.sim.dt_default = 1/30;
opt.nr.use_fsolve = true;
opt.verbose = true;

if LoadType
    ps.shunt(:,C.sh.frac_S)=0;
    ps.shunt(:,C.sh.frac_E)=1;
    fprintf('Load is an exponential type with lambda: 0.08\n')
    keyboard
else
    ps.shunt(:,C.sh.frac_S)=0;
    ps.shunt(:,C.sh.frac_E)=1;
    fprintf('Load is all constant P and Q \n')
    keyboard
end

% initialize the case
ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
% build the machine variables
[ps.mac,ps.exc,ps.gov] 		= get_mac_state(ps,'salient');
% initialize relays
ps.relay                    = get_relays(ps,'temperature',opt);

if check_derivatives    
    % build x and y
    [x0,y0] = get_xy(ps);
    % check df_dx
    x = x0; y = y0;
    disp('Checking df_dx for x0');
    [f0,df_dx,df_dy] = differential_eqs(0,x,y,ps,opt);
    f_handle = @(x) differential_eqs(0,x,y,ps,opt);
    checkDerivatives(f_handle,df_dx,x);
    % repeat from a perturbed starting point
    disp('Checking df_dx for perturbed x0');
    x = x0+randn(size(x0))*0.1;
    [f0,df_dx,df_dy] = differential_eqs(0,x,y,ps,opt);
    checkDerivatives(f_handle,df_dx,x);
    
    % check df_dy
    x = x0; y = y0;
    disp('Checking df_dy for y0');
    [f0,df_dx,df_dy] = differential_eqs(0,x,y,ps,opt);
    f_handle = @(y) differential_eqs(0,x,y,ps,opt);
    checkDerivatives(f_handle,df_dy,y);
    % repeat from a perturbed starting point
    disp('Checking df_dy for perturbed y0');
    y = y0+randn(size(y0))*0.1;
    [f0,df_dx,df_dy] = differential_eqs(0,x,y,ps,opt);
    checkDerivatives(f_handle,df_dy,y);
    
    % check dg_dx
    x = x0; y = y0;
    disp('Checking dg_dx for x0');
    [f0,dg_dx,dg_dy] = algebraic_eqs(0,x,y,ps,opt);
    g_handle = @(x) algebraic_eqs(0,x,y,ps,opt);
    checkDerivatives(g_handle,dg_dx,x);
    % repeat from a perturbed starting point
    disp('Checking dg_dx for perturbed x0');
    x = x0+randn(size(x0))*0.1;
    [f0,dg_dx,dg_dy] = algebraic_eqs(0,x,y,ps,opt);
    checkDerivatives(g_handle,dg_dx,x);
    
    % check dg_dy
    x = x0; y = y0;
    disp('Checking dg_dy for y0');
    [f0,dg_dx,dg_dy] = algebraic_eqs(0,x,y,ps,opt);
    g_handle = @(y) algebraic_eqs(0,x,y,ps,opt);
    checkDerivatives(g_handle,dg_dy,y);
    % repeat from a perturbed starting point
    disp('Checking dg_dy for perturbed y0');
    y = y0+randn(size(y0))*0.1;
    [f0,dg_dx,dg_dy] = algebraic_eqs(0,x,y,ps,opt);
    checkDerivatives(g_handle,dg_dy,y);
    keyboard
end

%% build an event matrix
event = zeros(4,C.ev.cols);
% start
event(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
% trip a branch
event(2,[C.ev.time C.ev.type]) = [10 C.ev.trip_branch];
event(2,C.ev.branch_loc) = 2;   % the open of branch 2, 3, or 6 will stop to integrate when it almost reaches the steady state
% set the end time
event(3,[C.ev.time C.ev.type]) = [t_max C.ev.finish];

%% run the simulation
[outputs,ps] = simgrid(ps,event,'sim_case9',opt);


%% print the results
fname = outputs.outfilename;
[t,delta,omega,Pm,Eap,temp,Vmag,theta,E1,Efd] = read_outfile(fname,ps);
omega_0 = 2*pi*ps.frequency;
omega_pu = omega / omega_0;

figure(1); clf; hold on; 
plot(t',omega_pu');
ylabel('\omega (pu)','FontSize',18);
xlabel('time (sec.)','FontSize',18);

figure(2); clf; hold on; 
plot(t',theta');
ylabel('\theta','FontSize',18);
xlabel('time (sec.)','FontSize',18);

figure(3); clf; hold on; 
plot(t',Vmag');
ylabel('|V|','FontSize',18);
xlabel('time (sec.)','FontSize',18);

figure(4); clf; hold on; 
plot(t',temp');
ylabel('temperature','FontSize',18);
xlabel('time (sec.)','FontSize',18);

figure(5); clf; hold on; 
plot(t',Pm');
ylabel('Pm','FontSize',18);
xlabel('time (sec.)','FontSize',18);

figure(6); clf; hold on; 
plot(t',delta'.*180./pi);
ylabel('Delta','FontSize',18);
xlabel('time (sec.)','FontSize',18);

figure(7); clf; hold on; 
plot(t',Eap');
ylabel('Eap','FontSize',18);
xlabel('time (sec.)','FontSize',18);

figure(8); clf; hold on; 
plot(t',E1');
ylabel('E1','FontSize',18);
xlabel('time (sec.)','FontSize',18);

figure(9); clf; hold on; 
plot(t',Efd');
ylabel('Efd','FontSize',18);
xlabel('time (sec.)','FontSize',18);
