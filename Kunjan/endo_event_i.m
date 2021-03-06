function [value,isterminal,direction] = endo_event_i(~,xy,~,ix,ps)
% usage: [value,isterminal,direction] = endo_event_i(~,xy,~,ix,ps)
% defines an exogenous event (relay) that will stop the DAE integration

% constants and settings
C = psconstants;
Vmag_threshold          = ps.relay(ix.re.uvls,C.re.threshold); % undervoltage threshold
omega_pu_threshold      = ps.relay(ix.re.ufls,C.re.threshold); % underfrequency threshold
dist_zone1_threshold    = ps.relay(ix.re.dist,C.re.threshold); % distance relay zone 1 setting

% extract some info from the inputs
x               = xy(1:ix.nx); 
y               = xy(ix.nx+1:end);
Vmags           = y(ix.y.Vmag);
Thetas          = y(ix.y.theta);
V               = Vmags .* exp(1i*Thetas);
If              = ps.Yf * V;
Imag_f          = abs(If);
temp_threshold  = ps.relay(ix.re.temp,C.re.threshold);
F               = ps.bus_i(ps.branch(:,C.br.from));
y_apparent      = Imag_f./Vmags(F);

% get the shunt_bus locations
sh_bus_ix = ps.bus_i(ps.shunt(:,1));
Vmag_sh   = y(ix.y.Vmag(sh_bus_ix));

% build zero crossing function
value = [temp_threshold       - x(ix.x.temp);       % zero crossing at maximum temperature
         Vmag_sh              - Vmag_threshold;     % trigger undervoltage load shedding
         x(ix.x.omega_pu)     - omega_pu_threshold; % trigger underfrequency load shedding
         dist_zone1_threshold - y_apparent];        % trigger zone 1 distance relay
     
% for relays that have already tripped set value to be in a safe range
is_tripped = ps.relay(:,C.re.tripped)==1;
value(is_tripped) = 10;
%value(ix.re.ufls) = 10;
%value(ix.re.dist) = 10;

% set directions
n_relays = size(ps.relay,1);
isterminal  =    ones(n_relays,1);    % terminate integration
direction   =   -ones(n_relays,1);    % to detect crossing from + to -
