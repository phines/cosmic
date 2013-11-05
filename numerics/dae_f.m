function [f,df_dx,df_dy] = dae_f(t,x,y)

% function
f = -y;

% derivative
df_dx = 0;
df_dy = -1;
