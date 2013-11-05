%% driver that tests the f and g equations for SMIB case
clear all; close all; clc; C = psconstants;

% initial conditions
Pm  = 1;
V2  = 1;
X12 = 0.1;
Xd  = 0.1;
Ea  = 1.1;
delta = asin(Pm*(X12+Xd)/Ea/V2);

x0 = [delta  0  Pm   Ea]';   % [delta delta_dot Pm Ea]
y0 = [1 0]';                 % [Vmag1 Theta1]


%% check derivatives for differential equations
%% check df_dx
x = x0; y = y0;
disp('Checking df_dx for x0');
[f0,df_dx,df_dy] = dae_smib_f(0,x,y);
f_handle = @(x) dae_smib_f(0,x,y);
checkDerivatives(f_handle,df_dx,x);
% repeat from a perturbed starting point
disp('Checking df_dx for perturbed x0');
x = x0+randn(size(x0))*0.1;
[f0,df_dx,df_dy] = dae_smib_f(0,x,y);
checkDerivatives(f_handle,df_dx,x);

%% check df_dy
x = x0; y = y0;
disp('Checking df_dy for y0');
[f0,df_dx,df_dy] = dae_smib_f(0,x,y);
f_handle = @(y) dae_smib_f(0,x,y);
checkDerivatives(f_handle,df_dy,y);
% repeat from a perturbed starting point
disp('Checking df_dy for perturbed y0');
y = y0+randn(size(y0))*0.1;
[f0,df_dx,df_dy] = dae_smib_f(0,x,y);
checkDerivatives(f_handle,df_dy,y);

%% check dg_dx
x = x0; y = y0;
disp('Checking dg_dx for x0');
[f0,dg_dx,dg_dy] = dae_smib_g(0,x,y);
g_handle = @(x) dae_smib_g(0,x,y);
checkDerivatives(g_handle,dg_dx,x);
% repeat from a perturbed starting point
disp('Checking dg_dx for perturbed x0');
x = x0+randn(size(x0))*0.1;
[f0,dg_dx,dg_dy] = dae_smib_g(0,x,y);
checkDerivatives(g_handle,dg_dx,x);

%% check dg_dy
x = x0; y = y0;
disp('Checking dg_dy for y0');
[f0,dg_dx,dg_dy] = dae_smib_g(0,x,y);
g_handle = @(y) dae_smib_g(0,x,y);
checkDerivatives(g_handle,dg_dy,y);
% repeat from a perturbed starting point
disp('Checking dg_dy for perturbed y0');
y = y0+randn(size(y0))*0.1;
[f0,dg_dx,dg_dy] = dae_smib_g(0,x,y);
checkDerivatives(g_handle,dg_dy,y);

return
