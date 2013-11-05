function [f,df_dx,df_dy] = dae_smib_f(t,x,y)
% test differential equation for a single machine infinite bus

% constants
D       = 1.5;
M       = 3;
Xd      = 0.1;
ix.nx   = 4;
ix.ny   = 2;

% extract variables
delta       = x(1);
delta_dot   = x(2);
Pm          = x(3);
Ea          = x(4);    
Vmag1       = y(1);
Theta1      = y(2);

% initialize the output
f = zeros(length(x),1);

% calculate Pe
Pe = ((Ea.*Vmag1)./Xd).*sin(delta-Theta1); 

% calculate swing equations
f(1) = delta_dot;
f(2) = (Pm - Pe - D.*delta_dot)./M;
%f(3) = 0.008;

% output df_dx and df_dy if requested
if nargout>1

    % build df_dx
    dPg_ddelta = (Ea.*Vmag1./Xd).*cos(delta-Theta1);
    
    dFswing_ddelta_values = -dPg_ddelta./M;
    dFswing_domega_values = -D./M;
    dFswing_dPm_values    = 1./M;
    dFswing_dEa_values    = -sin(delta-Theta1).*Vmag1./(M.*Xd);

    % assemble df_dx
    df_dx = sparse(ix.nx,ix.nx);
    
    % dFdelta_dot_domega
    df_dx = df_dx + sparse(1,2,1,ix.nx,ix.nx);
    
    % dFswing_ddelta
    df_dx = df_dx + sparse(2,1,dFswing_ddelta_values,ix.nx,ix.nx);
    % dFswing_domega
    df_dx = df_dx + sparse(2,2,dFswing_domega_values,ix.nx,ix.nx);
    % dFswing_dPm
    df_dx = df_dx + sparse(2,3,dFswing_dPm_values,ix.nx,ix.nx);
    % dFswing_dEa
    df_dx = df_dx + sparse(2,4,dFswing_dEa_values,ix.nx,ix.nx);
    
end
if nargout>2
    
    df_dy = sparse(ix.nx,ix.ny);
 
    % build df_dy (change in f wrt the algebraic variables)
    dPg_dtheta = Ea.*sin(delta-Theta1)./Xd;
    dFswing_dVmag_values = -(dPg_dtheta)./M;
    
	% assemble df_dy
    df_dy = df_dy + sparse(2,1,dFswing_dVmag_values,ix.nx,ix.ny);
    df_dy = df_dy + sparse(2,2,-dFswing_ddelta_values,ix.nx,ix.ny);
    
end
