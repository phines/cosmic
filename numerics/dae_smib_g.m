function [g,dg_dx,dg_dy] = dae_smib_g(t,x,y)
% test algebraic equation for a single machine infinite bus

% constants
Vmag2   = 1;
X12     = 0.1;
Xd      = 0.1;
ix.nx   = 4;
ix.ny   = 2;

% extract variables
delta  = x(1);
Ea     = x(4);   
Vmag1  = y(1);
Theta1 = y(2);

% initialize output
g = zeros(length(y),1);

% calculate power from bus 1 to bus 2 
P12     = +(Vmag1*Vmag2/X12)*sin(Theta1); % Paul changed this to plus
Q12     = -(Vmag1*Vmag2/X12)*cos(Theta1) + (Vmag1^2)/X12; % Paul changed signs here

% calculate P and Q for a classic generator model
P1g = +(Ea*Vmag1/Xd)*sin(Theta1-delta);% Paul changed signs here
Q1g = -(Ea*Vmag1/Xd)*cos(Theta1-delta) + (Vmag1^2)/Xd; % Paul changed signs here

% output
g(1) = P12 + P1g;
g(2) = Q12 + Q1g;

% output dg_dx and dg_dy if requested
if nargout>1
    % active power
    dP1g_ddelta_values = -(Ea.*Vmag1./Xd) .* cos(Theta1-delta);
    dP1g_dEa_values    = (Vmag1./Xd) .* sin(Theta1-delta);
    
    % reactive power
    dQ1g_ddelta_values = +(Ea.*Vmag1./Xd) .* (-sin(Theta1-delta));
    dQ1g_dEa_values    = -(Vmag1./Xd) .* cos(Theta1-delta);
    
    % assemble dg_dx
    dg_dx = sparse(ix.ny,ix.nx);
    % dPg_ddelta
    dg_dx = dg_dx + sparse(1,1,dP1g_ddelta_values,ix.ny,ix.nx);
    % dPg_dEa
    dg_dx = dg_dx + sparse(1,4,dP1g_dEa_values,ix.ny,ix.nx);
    % dQg_ddelta
    dg_dx = dg_dx + sparse(2,1,dQ1g_ddelta_values,ix.ny,ix.nx);
    % dQg_dEa
    dg_dx = dg_dx + sparse(2,4,dQ1g_dEa_values,ix.ny,ix.nx);
end
if nargout>2
    % assemble dg_dy
    dg_dy = sparse(ix.ny,ix.ny);
    %P12     = +(Vmag1*Vmag2/X12)*sin(Theta1); % Paul changed this to plus
    dP12_dVmag1  = +Vmag2/X12*sin(Theta1);
    dP12_dTheta1 = (Vmag1*Vmag2/X12)*cos(Theta1);
    
    %Q12     = -(Vmag1*Vmag2/X12)*cos(Theta1) + (Vmag1^2)/X12; % Paul changed signs here
    dQ12_dVmag1  = -(Vmag2/X12)*cos(Theta1) + (2*Vmag1)/X12;
    dQ12_dTheta1 = +(Vmag1*Vmag2/X12)*sin(Theta1);
    
    % assemble dP_IK_dVmag_I and K
    dg_dy = dg_dy + sparse(1,1,dP12_dVmag1,ix.ny,ix.ny);
    % assemble dQ_IK_dVmag_I and K
    dg_dy = dg_dy + sparse(2,1,dQ12_dVmag1,ix.ny,ix.ny);
    % assemble dP_IK_dtheta_I and K
    dg_dy = dg_dy + sparse(1,2,dP12_dTheta1,ix.ny,ix.ny);
    % assemble dQ_IK_dtheta_I and K
    dg_dy = dg_dy + sparse(2,2,dQ12_dTheta1,ix.ny,ix.ny);
    
    % derivatives for real power injection from generators
    %P1g = +(Ea*Vmag1/Xd)*sin(Theta1-delta);
    dP1g_dVgen = +(Ea/Xd)*sin(Theta1-delta);
    dP1g_dtheta = +(Ea*Vmag1/Xd)*cos(Theta1-delta);
    % derivatives for reactive power injection from generators
    %Q1g = -(Ea*Vmag1/Xd)*cos(Theta1-delta) + (Vmag1^2)/Xd;
    dQ1g_dVgen = -(Ea/Xd)*cos(Theta1-delta) + (2*Vmag1)/Xd;
    dQ1g_dtheta = -(Ea*Vmag1/Xd)*(-sin(Theta1-delta));
    
    % stick the stuff above into dg_dy
    dg_dy = dg_dy + sparse(1,1,dP1g_dVgen,ix.ny,ix.ny);
    dg_dy = dg_dy + sparse(2,1,dQ1g_dVgen,ix.ny,ix.ny);
    dg_dy = dg_dy + sparse(1,2,dP1g_dtheta,ix.ny,ix.ny);
    dg_dy = dg_dy + sparse(2,2,dQ1g_dtheta,ix.ny,ix.ny);
end
