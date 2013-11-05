% test code for dae_trap
clear all; 

% time span 
t_max = 20;
t_span = [0 t_max];

% initial conditions
x0 = sin(0);
y0 = -cos(0);

% set up the differential function
f = @(t,x,y) dae_f(t,x,y);
% set up the algebraic function
g = @(t,x,y) dae_g(t,x,y);

% integrate
[t_steps,X,Y] = solve_dae(f,g,x0,y0,t_span); 

% plot X and Y
figure(1); clf; 
plot(t_steps,X,t_steps,Y); hold on;
plot(t_steps,-cos(t_steps),'r');
plot(t_steps,sin(t_steps),'c');
legend('x(t)','y(t)','x_true','y_true');
xlabel('time');
axis([0 t_max -inf inf]);

% plot time steps
figure(2); clf;
dt = t_steps(2:end) - t_steps(1:(end-1));
plot(t_steps(1:(end-1)),dt);
xlabel('time');
axis([0 t_max -inf inf]);
%set(gca,'yscale','log')
