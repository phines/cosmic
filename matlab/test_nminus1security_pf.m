%%  test N-1 security in terms of ac power flow
clear all; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end
outputs = []; psout = [];

load case2383_mod_ps_dyn;
ps.branch(:,C.br.tap)       = 1;
ps.shunt(:,C.sh.frac_S)     = 0;
ps.shunt(:,C.sh.frac_E)     = 1;
ps.shunt(:,C.sh.frac_Z)     = 0;
ps.shunt(:,C.sh.gamma)      = 0.08;

% ps = updateps(case39_ps);
ps = newpf(ps);
opt = psoptions;
opt.verbose = false;
n_br  = size(ps.branch,1);
for i = 1:n_br
    fprintf('running contingency #%d: branches %d\n',i,i);
    ps.branch(i,C.br.status) = 0;  % contingency no. i
    [ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
    % build the machine variables
    [ps.mac,ps.exc,ps.gov] 		= get_mac_state(ps,'salient');
    % build indices
    n  = size(ps.bus,1);
    ng = size(ps.gen,1);
    m  = size(ps.branch,1);
    n_sh = size(ps.shunt,1);
    ix   = get_indices(n,ng,m,n_sh,opt);
    % build x and y
    [x0,y0] = get_xy(ps,opt);
    % solve algebraic only
    y_new = solve_algebraic([],x0,y0,ps,opt);
    
    if isempty(y_new)
        psout{i}= [];
        outputs{i} = false;
    else
        psout{i}= y_new;
        outputs{i} = true;
    end
        
    ps.branch(i,C.br.status) = 1;  % put no.i branch back to service    
end
