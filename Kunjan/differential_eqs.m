function [f,df_dx,df_dy] = differential_eqs(t,x,y,ps,opt) %#ok<INUSL>
% usage: [f,df_dx,df_dy] = differential_eqs(t,x,y,ps,opt)
% differential equations that model the elements of the power system
%
% inputs:
%   t   -> time in seconds
%   x   -> [delta omega Pm Eap E1 Efd] for each machine and [temp] for each branch
%   y   -> [Vmag Theta]
%   ps  -> power system structure
%   opt -> options structure
%
% outputs: 
%   f(1) -> ddelta_dt
%   f(2) -> domega_dt
%   f(3) -> dPm_dt
%   f(4) -> dEap_dt
%   f(5) -> dE1_dt
%   f(6) -> dEfd_dt
%   f(7) -> dtemp_dt

% constants
C           = psconstants;
n           = size(ps.bus,1);
ng          = size(ps.mac,1);
m           = size(ps.branch,1);
n_sh        = size(ps.shunt,1);
ix          = get_indices(n,ng,m,n_sh);

% extract parameters from ps
Xds     = ps.mac(:,C.ma.Xd);
Xdps    = ps.mac(:,C.ma.Xdp);
Xqs     = ps.mac(:,C.ma.Xq);
Td0ps   = ps.mac(:,C.ma.Td0p);
Ds      = ps.mac(:,C.ma.D);
Ms      = ps.mac(:,C.ma.M);
omega_0 = 2*pi*ps.frequency;

% extract differential variables
deltas      = x(ix.x.delta);
omegas_pu   = x(ix.x.omega_pu);
Pms         = x(ix.x.Pm);
Eaps        = x(ix.x.Eap);
Efds        = x(ix.x.Efd);
E1s         = x(ix.x.E1);
temps       = x(ix.x.temp);

% extract algebraic variables and fix the slack bus angle
mac_buses   = ps.mac(:,C.mac.gen);
mac_bus_i   = ps.bus_i(mac_buses);

Vmags       = y(ix.y.Vmag);
Thetas      = y(ix.y.theta);
mac_Vmags   = Vmags(mac_bus_i);
mac_Thetas  = Thetas(mac_bus_i);

% machine angles, relative to the bus angles
delta_sys   = y(ix.y.delta_sys);
delta_m     = deltas + delta_sys - mac_Thetas;

% calculate Pe: Eq. 7.81 from Bergen & Vittal
Pes = (Eaps.*mac_Vmags./Xdps).*sin(delta_m) + mac_Vmags.^2./2.*(1./Xqs-1./Xdps).*sin(2*delta_m);

% initialize output
f = zeros(length(x),1);

% calculate swing equations
f(ix.f.delta_dot) = omega_0.*(omegas_pu-1);
f(ix.f.omega_dot) = (Pms - Pes - Ds.*(omegas_pu-1))./Ms;
f(ix.f.Eap_dot)   = -Eaps.*Xds./(Td0ps.*Xdps)+(Xds./Xdps-1).*mac_Vmags.*cos(delta_m)./Td0ps+Efds./Td0ps; % Eq. 7.75 from Bergen & Vittal   
[f(ix.f.Pm_dot),df_dx_gov]                              = governor_eqs(Pms',omegas_pu,ps);
[f(ix.f.Efd_dot),f(ix.f.E1_dot),df_dx_exc,df_dy_exc]   	= exciter_eqs([Efds';E1s'],mac_Vmags,ps);

% calculate the temperature change
[f(ix.f.temp_dot),df_dy_temp]   = temp_eqs(Vmags,Thetas,temps,ps,opt);

% output df_dx and df_dy if requested
if nargout>1
    % build df_dx
    dPg_ddelta = (Eaps.*mac_Vmags./Xdps).*cos(delta_m) + mac_Vmags.^2.*(1./Xqs-1./Xdps).*cos(2*delta_m);
    dFswing_ddelta_values = -dPg_ddelta./Ms;
    dFswing_domega_values = -Ds./Ms;
    dFswing_dPm_values    = 1./Ms;
    dFswing_dEa_values    = -sin(delta_m).*mac_Vmags./(Ms.*Xdps);

    dPm_domegas_values     = df_dx_gov(:,1);
    dPm_dPm_values         = df_dx_gov(:,2);
    dEap_dEap_values       = -Xds./(Td0ps.*Xdps);
    dEap_ddelta_values     = -(Xds./Xdps-1).*mac_Vmags.*sin(delta_m)./Td0ps;
    dEap_dEfd_values       = 1./Td0ps;
    dE1_dE1_values         = df_dx_exc(:,1);
    dEfd_dE1_values        = df_dx_exc(:,2);
    dEfd_dEfd_values       = df_dx_exc(:,3);    

    % assemble df_dx
    df_dx = sparse(ix.nx,ix.nx);
    % dFswing_ddelta
    df_dx = df_dx + sparse(ix.f.omega_dot,ix.x.delta,dFswing_ddelta_values,ix.nx,ix.nx);
    % dFswing_domega
    df_dx = df_dx + sparse(ix.f.omega_dot,ix.x.omega_pu,dFswing_domega_values,ix.nx,ix.nx);
    % dFswing_dPm
    df_dx = df_dx + sparse(ix.f.omega_dot,ix.x.Pm,dFswing_dPm_values,ix.nx,ix.nx);
    % dFswing_dEa
    df_dx = df_dx + sparse(ix.f.omega_dot,ix.x.Eap,dFswing_dEa_values,ix.nx,ix.nx);
    % dFdelta_dot_domega
    df_dx = df_dx + sparse(ix.f.delta_dot,ix.x.omega_pu,omega_0,ix.nx,ix.nx);
    % dFtemp_dot_dtemp
    df_dx = df_dx + sparse(ix.f.temp_dot,ix.x.temp,-opt.sim.temp.K,ix.nx,ix.nx);
	
	% dPm_domega
    df_dx = df_dx + sparse(ix.f.Pm_dot,ix.x.omega_pu,dPm_domegas_values,ix.nx,ix.nx);
    % dPm_dot_dPm
    df_dx = df_dx + sparse(ix.f.Pm_dot,ix.x.Pm,dPm_dPm_values,ix.nx,ix.nx);
    % dEap_dot_dEap
    df_dx = df_dx + sparse(ix.f.Eap_dot,ix.x.Eap,dEap_dEap_values,ix.nx,ix.nx);
    % dEap_dot_ddelta
    df_dx = df_dx + sparse(ix.f.Eap_dot,ix.x.delta,dEap_ddelta_values,ix.nx,ix.nx);
    % dEap_dot_dEfd
    df_dx = df_dx + sparse(ix.f.Eap_dot,ix.x.Efd,dEap_dEfd_values,ix.nx,ix.nx);
    % dE1_dot_dE1
    df_dx = df_dx + sparse(ix.f.E1_dot,ix.x.E1,dE1_dE1_values,ix.nx,ix.nx);
    % dEfd_dot_dE1
    df_dx = df_dx + sparse(ix.f.Efd_dot,ix.x.E1,dEfd_dE1_values,ix.nx,ix.nx);
    % dEfd_dot_dEfd
    df_dx = df_dx + sparse(ix.f.Efd_dot,ix.x.Efd,dEfd_dEfd_values,ix.nx,ix.nx);
end
if nargout>2
    % build df_dy (change in f wrt the algebraic variables)
    dPg_dtheta = Eaps.*sin(delta_m)./Xdps + mac_Vmags.*(1./Xqs-1./Xdps).*sin(2*delta_m);
    dFswing_dVmag_values        = -(dPg_dtheta)./Ms;
	dEap_dVmag_values           = (Xds./Xdps-1).*cos(delta_m)./Td0ps;
    dE1_dVmag_values            = df_dy_exc(:,1);
    dEfd_dVmag_values           = df_dy_exc(:,2);
    dFtemp_dot_dVmag_f_values   = df_dy_temp(:,1);
    dFtemp_dot_dVmag_t_values   = df_dy_temp(:,2);
    dFtemp_dot_dtheta_values    = df_dy_temp(:,3); 
    
	% assemble df_dy
    cols = ix.y.Vmag(mac_bus_i);
    df_dy = sparse(ix.f.omega_dot,cols,dFswing_dVmag_values,ix.nx,ix.ny);
    cols = ix.y.theta(mac_bus_i);
    df_dy = df_dy + sparse(ix.f.omega_dot,cols,-dFswing_ddelta_values,ix.nx,ix.ny);
    cols = ix.y.delta_sys;
    df_dy = df_dy + sparse(ix.f.omega_dot,cols,dFswing_ddelta_values,ix.nx,ix.ny);
	% dEap_dot_dVmag
    cols  = ix.y.Vmag(mac_bus_i);
    df_dy = df_dy + sparse(ix.f.Eap_dot,cols,dEap_dVmag_values,ix.nx,ix.ny);
    % dEap_dot_ddelta_sys
    cols = ix.y.delta_sys;
    df_dy = df_dy + sparse(ix.f.Eap_dot,cols,dEap_ddelta_values,ix.nx,ix.ny);
    % dEap_dot_dtheta
    cols = ix.y.theta(mac_bus_i);
    df_dy = df_dy + sparse(ix.f.Eap_dot,cols,-dEap_ddelta_values,ix.nx,ix.ny);
    % dE1_dot_dVmag
    cols  = ix.y.Vmag(mac_bus_i);
    df_dy = df_dy + sparse(ix.f.E1_dot,cols,dE1_dVmag_values,ix.nx,ix.ny);
    % dEfd_dot_dVmag
    cols  = ix.y.Vmag(mac_bus_i);
    df_dy = df_dy + sparse(ix.f.Efd_dot,cols,dEfd_dVmag_values,ix.nx,ix.ny);
    
    % dFtemp_dot
    br_status  = (ps.branch(:,C.br.status)>=1);
    F       = ps.bus_i(ps.branch(br_status,C.br.from));
    T       = ps.bus_i(ps.branch(br_status,C.br.to));
    
    cols  = ix.y.Vmag(F);   % insert the two components of the F side
    df_dy = df_dy + sparse(ix.f.temp_dot(br_status),cols,dFtemp_dot_dVmag_f_values(br_status),ix.nx,ix.ny);
    cols  = ix.y.theta(F);   
    df_dy = df_dy + sparse(ix.f.temp_dot(br_status),cols,dFtemp_dot_dtheta_values(br_status),ix.nx,ix.ny);
    
    cols  = ix.y.Vmag(T);  % insert the two components of the T side
    df_dy = df_dy + sparse(ix.f.temp_dot(br_status),cols,dFtemp_dot_dVmag_t_values(br_status),ix.nx,ix.ny);
    cols  = ix.y.theta(T);
    df_dy = df_dy + sparse(ix.f.temp_dot(br_status),cols,-dFtemp_dot_dtheta_values(br_status),ix.nx,ix.ny);
end
