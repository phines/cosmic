function [value,isterminal,direction] = endo_event_rec(t,xy,ix,ps,dt)
% usage: [value,isterminal,direction] = endo_event(~,xy,ix,ps)
% defines an exogenous event (relay) that will stop the DAE integration

% constants and settings
global dist2threshold state_a 
C = psconstants;
oc_setting1             = ps.relay(ix.re.oc,C.re.setting1);    % over-current maximum current
temp_setting1           = ps.relay(ix.re.temp,C.re.setting1);  % 
temp_K                  = ps.relay(ix.re.temp,C.re.setting2);
Vmag_threshold          = ps.relay(ix.re.uvls,C.re.threshold); % undervoltage threshold
omega_pu_threshold      = ps.relay(ix.re.ufls,C.re.threshold); % underfrequency threshold
dist_zone1_threshold    = ps.relay(ix.re.dist,C.re.threshold); % distance relay zone 1 setting
SMALL_EPS = 1e-12;

% extract some info from the inputs
x               = xy(1:ix.nx); 
y               = xy(ix.nx+1:end);
Vr              = y(ix.y.Vr);
Vi              = y(ix.y.Vi);
V               = Vr + 1i*Vi;
Vmags           = abs(V);
If              = ps.Yf * V;
It              = ps.Yt * V;
Imag_f          = abs(If);
Imag_t          = abs(It);
Imag            = max(Imag_f, Imag_t);
temp            = Imag.^2 ./ temp_K .* (1 - exp(-temp_K*t));

if ~isempty(dt)
    state_a(ps.relay(ix.re.temp,C.re.id)) = max(state_a(ps.relay(ix.re.temp,C.re.id))+(temp-temp_setting1.^2./temp_K)*dt+SMALL_EPS, 0);
    state_a(ps.relay(ix.re.oc,C.re.id)) = max(state_a(ps.relay(ix.re.oc,C.re.id))+(Imag-oc_setting1)*dt+SMALL_EPS, 0);
end

dist2threshold(ps.relay(ix.re.temp,C.re.id)) = ...
    ps.relay(ix.re.temp,C.re.threshold) - state_a(ps.relay(ix.re.temp,C.re.id)); % associate them with their global id
dist2threshold(ps.relay(ix.re.oc,C.re.id)) = ...
    ps.relay(ix.re.oc,C.re.threshold) - state_a(ps.relay(ix.re.oc,C.re.id)); % associate them with their global id
dist2threshold(dist2threshold<0)=0;

F               = ps.bus_i(ps.branch(:,C.br.from));
y_apparent      = Imag_f./Vmags(F);
nload           = size(ps.shunt,1);
load_freq       = zeros(nload,1);
near_gen        = ps.shunt(:,C.sh.near_gen);
for i = 1:nload
    near_gen_id = ps.gen_i(near_gen(i));
    if near_gen_id ~= 0
        load_freq_source = ix.x.omega_pu(near_gen_id);
        load_freq(i)  = x(load_freq_source);
    end
end

% get the shunt_bus locations
sh_bus_ix = ps.bus_i(ps.shunt(:,1));
Vmag_sh   = Vmags(sh_bus_ix);

% build zero crossing function
value = [temp_setting1.^2./temp_K - temp;           % active temperature relay
         oc_setting1          - Imag;               % active over current relay
         Vmag_sh              - Vmag_threshold;     % active undervoltage load shedding
         load_freq            - omega_pu_threshold; % active underfrequency load shedding
         dist_zone1_threshold - y_apparent];        % active zone 1 distance relay

% for relays that have already tripped set value to be in a safe range
is_tripped = ps.relay(:,C.re.tripped)==1;
value(is_tripped) = 10;
%value(ix.re.ufls) = 10;
%value(ix.re.dist) = 10;

% set directions
n_relays = size(ps.relay,1);
isterminal  =    ones(n_relays,1);    % terminate integration
direction   =   -ones(n_relays,1);    % to detect crossing from + to -
