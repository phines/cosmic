function [outputs] = get_event_thresholds(inputs)
% usage: [outputs] = get_event_thresholds(inputs)
%
% calculate distance to discrete changes for dae solver
% currently, these events come from power system relays

C               = psconstants;       % power system constants
relay_event     = [];
F               = ps.bus_i(ps.branch(:,C.br.from));

% relay info
dist_relays          = ps.relay(ix.re.dist,:);
dist_zone1_threshold = ps.relay(ix.re.dist,C.re.threshold); % distance relay zone 1 setting

Vmags1          = y1(ix.y.Vmag);
Thetas1         = y1(ix.y.theta);
V1              = Vmags1 .* exp(1i*Thetas1);
If1             = ps.Yf * V1;
Imag_f1         = abs(If1);
y_apparent1     = Imag_f1./Vmags1(F);
relay_mismatch1 = dist_zone1_threshold - y_apparent1
if any(relay_mismatch1 < 0) % zero crossing, + to -
    fprintf('  t = %.2f: Distance relay at branch %d tripped...\n',t1,dist_relays(relay_mismatch1<0,C.re.branch_loc));
    keyboard
end

% fprintf('Algorithm failure in the newton step, assume voltage collapse.\n');
Vmags       = Y(ix.y.Vmag,end);
relay_event = Vmags < (median(Vmags) - 2.*std(Vmags));