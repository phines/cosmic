function [mac,exc,gov] = get_mac_state(ps,mode)
% usage: [mac,exc,gov] = get_mac_state(ps,mode)
%  ps - ps structure typical
%  mode is one of the following:
%   'classical', Ea behind Xd
%   'salient', Salient pole machine

% NOTE that this function assumes that there is one mac entry per gen, in
% the same sequence.

% grab some data
C  = psconstants;
j  = 1i;
Pg = ps.gen(:,C.ge.Pg)/ps.baseMVA;
Qg = ps.gen(:,C.ge.Qg)/ps.baseMVA;
G  = ps.bus_i(ps.gen(:,1));
theta_g     = ps.bus(G,C.bu.Vang)*pi/180;
Va          = ps.bus(G,C.bu.Vmag).*exp(1i*theta_g);
Sa          = Pg + j*Qg; 
Ia          = conj(Sa./Va);
phi         = angle(Sa);        % power factor angle for the generator

if abs(phi) > pi/2
    error('strange angle for the power factor');
end
Xd  = ps.mac(:,C.ma.Xd);
Xdp = ps.mac(:,C.ma.Xdp);
Xq  = ps.mac(:,C.ma.Xq);
mac = ps.mac;
exc = ps.exc;
gov = ps.gov;
ng = size(mac,1);
omega_0 = 2*pi*ps.frequency;

switch mode
    case 'classical'
        % find Pm in per unit
        mac(:,C.ma.Pm) = Pg; % at nominal frequency
        % find delta
        delta_m_theta = atan2( Pg , ( Qg + Va.^2 ./ Xd ) );
        delta = delta_m_theta; % changed from delta_m_theta + theta_g;
        mac(:,C.ma.delta) = delta;
        % find Ea
        Ea_mag = Pg .* Xd ./ ( Va .* sin(delta) );
        mac(:,C.ma.Ea) = Ea_mag;
        % delta_dot = 0
        mac(:,C.ma.omega) = omega_0;
    case 'salient'
        if any(Xq==0)
            error('we need Xq values to do salient pole model');
        end
        for i=1:ng
            % find aprime
            aprime = Va(i) + j*Xq(i)*Ia(i); % + r*Ia
            delta  = angle(aprime);
            delta_m = delta - theta_g(i);
            Vaq = abs(Va(i))*cos(delta_m);
            Iad_mag = abs(Ia(i)).*sin(delta-angle(Ia(i)));
            Eap_mag = Vaq + Xdp(i)*Iad_mag; % + r*Iad
            Ea_mag  = (Xd(i)-Xdp(i))*Iad_mag + Eap_mag;
            % verify the machine
            if ~verify_mac(delta_m, Ea_mag, Eap_mag, ...
                    Xd(i), Xdp(i), Xq(i), abs(Va(i)), Pg(i), Qg(i))
                error('machine %d salient model is not consistent...\n',i);
            end
            % save the results
            mac(i,C.ma.Ea)      = Ea_mag;
            mac(i,C.ma.Eap)     = Eap_mag;
            mac(i,C.ma.delta_m) = delta_m;
            mac(i,C.ma.omega)   = omega_0;
            mac(i,C.ma.Pm0)     = Pg(i); % at nominal frequency
            exc(i,C.ex.Vref)    = Ea_mag/ ps.exc(i,C.ex.Ke) + abs(Va(i));    
            exc(i,C.ex.E1)      = Ea_mag/ ps.exc(i,C.ex.Ke);                 
            exc(i,C.ex.Efd)     = Ea_mag;                                    
            gov(i,C.go.Pref)    = Pg(i);                                     
        end
    otherwise
        error('not a valid mode for get_mac_state');
end
