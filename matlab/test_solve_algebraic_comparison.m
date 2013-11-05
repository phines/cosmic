
clear all; close all; clc; C = psconstants;

if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end

%% select data case to simulate
% ps = updateps(case3_ps);           
ps = updateps(case9_ps);                   
% ps = updateps(case39_ps);           

% ps = replicate_case(ps,4);          

ps = unify_generators(ps);            % converts multiple generators connected to a bus into a single generator

%% solve power flow in polar form and rectangular form

% set options
opt = psoptions;
% opt.nr.use_fsolve = true;
% opt.pf.linesearch = 'cubic_spline';

[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);

% build the machine variables
[ps.mac,ps.exc,ps.gov] 		= get_mac_state(ps,'salient');   
% initialize relays
ps.relay                    = get_relays(ps,'temperature',opt);

t = 1;
C           = psconstants;
n           = size(ps.bus,1);
n_macs      = size(ps.mac,1);
m           = size(ps.branch,1);
n_shunts    = size(ps.shunt,1);
j           = 1i;

fprintf('Initialize ps \n')
% V       = ps.bus(:,C.bu.Vmag).*exp(1i * ps.bus(:,C.bu.Vang));
% ps.bus(:,C.bu.Vr) = real(V);
% ps.bus(:,C.bu.Vi) = imag(V);
ps = newpf_rec(ps,opt);

fprintf('Compare solve_algebraic \n')
polar = 0;
rec   = 1;

if rec
    ix_rec      = get_indices_rec(n,n_macs,m,n_shunts,opt);
    [x_rec,y_rec] = get_xy_rec(ps,opt);
    y_new_rec = solve_algebraic_rec(t,x_rec,y_rec,ps,opt);
%     keyboard
    Vr = y_new_rec(ix_rec.y.Vr);
    Vi = y_new_rec(ix_rec.y.Vi);
    Vmag_rec = abs(Vr + j.*Vi);
    Vang_rec = angle(Vr + j.*Vi);
    Vmag_Vang_rec = [Vmag_rec,Vang_rec];
end

if polar
    ix_polar    = get_indices(n,n_macs,m,n_shunts,opt);
    [x_polar,y_polar] = get_xy(ps,opt);
    y_new_polar = solve_algebraic(t,x_polar,y_polar,ps,opt);
    Vmag_polar = y_new_polar(ix_polar.y.Vmag);
    Vang_polar = y_new_polar(ix_polar.y.theta);
    Vmag_Vang_polar = [Vmag_polar,Vang_polar];
end


if polar && rec
    % compare the results from the above two functions
    err_Vmag_Vang = Vmag_Vang_rec - Vmag_Vang_polar;
    MaxErrorVmagVang = max(abs(err_Vmag_Vang));
    fprintf('The maximum discrepency for the Vmag and Vang results are %d and %d \n',MaxErrorVmagVang);
end

