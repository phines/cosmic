% test code for dae_smib_trap
clear all; close all; clc; 

% options
options = numerics_options;

% time span 
t_max = 20;
t1  = 5; % exogenous event
t_span1 = [0 t1];
t_span2 = [t1+options.sim.eps_thresh t_max];

% initial conditions
Pm  = 1;
V2  = 1;
X12 = 0.1;
Xd  = 0.1;
Ea  = 1.1;
delta = asin(Pm*(X12+Xd)/Ea/V2);

x0 = [delta  0  Pm   Ea]';   % [delta delta_dot Pm Ea]
y0 = [1 0]';                 % [Vmag1 Theta1]

% set up the differential function
f = @(t,x,y) dae_smib_f(t,x,y);
% set up the algebraic function
g = @(t,x,y) dae_smib_g(t,x,y);
% set up the event function
h = [];

% integrate first interval

[t_steps1,X1,Y1] = solve_dae(f,g,h,x0,y0,t_span1,options); 

x1 = X1(:,end); 
y1 = Y1(:,end);

x1(3) = 2; % increase Pm

[t_steps2,X2,Y2] = solve_dae(f,g,h,x1,y1,t_span2,options); 

figure(1); clf; 
plot([t_steps1;t_steps2],[X1 X2],[t_steps1;t_steps2],[Y1 Y2]); legend('x(t)','y(t)'); hold on;
xlabel('time');
axis([0 t_max -inf inf]);
return