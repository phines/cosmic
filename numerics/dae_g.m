function [g,dg_dx,dg_dy] = dae_g(t,x,y)

% function
g = sin(t) - x;

% derivative
dg_dx = -1;
dg_dy = 0;
