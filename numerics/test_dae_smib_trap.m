% test code for dae_smib_trap
clear all; close all; clc; 

% time span 
t_max = 19.345728;
t_span = [0 t_max];

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
h = @(t,x,y) dae_test_h(t,x,y);

% integrate
options = numerics_options;
[t_steps,X,Y] = solve_dae(f,g,h,x0,y0,t_span,options); 

figure(1); clf; 
plot(t_steps,X,t_steps,Y); legend('x(t)','y(t)'); hold on;
%plot(t_steps,-cos(t_steps),'r');
xlabel('time');
axis([0 t_max -inf inf]);
return
figure(2); clf;
dt = t_steps(2:end) - t_steps(1:(end-1));
plot(t_steps(1:(end-1)),dt);
xlabel('time');
axis([0 t_max -inf inf]);
%set(gca,'yscale','log')