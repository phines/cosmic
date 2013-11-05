
% Test newRpf funtion, and compare the results with the ones form newpf.m

clear all; close all; clc; C = psconstants;

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

if ~(ismcc || isdeployed)
    addpath('../data');
end

%% select data case to simulate
% ps = updateps(case3_ps);            % GOOD
% ps = updateps(case6_ps);            % GOOD
% ps = updateps(case9_ps);            % GOOD 
% ps = updateps(case14_ps);           % GOOD
% ps = updateps(case24_ps);           % GOOD --- This case has multiple generators at each gen bus.
% ps = updateps(case30_ps);           % GOOD
% ps = updateps(case39_ps);           % GOOD
% ps = updateps(case300_ps);          % GOOD
% 
% data = load('case2383_mod_ps.mat'); % GOOD
% ps = data.ps;
% 
data = load('case3120sp_ps.mat');   % GOOD
ps = data.ps;

% ps = replicate_case(ps,6);          % GOOD --- replicate case9_ps to form a large case 

ps = unify_generators(ps);            % converts multiple generators connected to a bus into a single generator

ps.branch(3,C.br.status) = 0;         % randomly picked it up
%% solve power flow in polar form and rectangular form

% set options
opt = psoptions;
% opt.pf.linesearch = 'cubic_spline';

polar = 1;
rec   = 1;

if polar
    tic
    ps_polar = newpf(ps,opt);  
    toc
end


if rec
    tic
    ps_rec = newpf_rec(ps,opt);
    toc
end


if polar && rec
    % compare the results from the above two functions
    err_Vmag_Vang = ps_polar.bus(:,8:9)-ps_rec.bus(:,8:9);
    err_Pg_Qg = ps_polar.gen(:,2:3)-ps_rec.gen(:,2:3);
    
    MaxErrorVmagVang = max(abs(err_Vmag_Vang));
    MaxErrorPgQg = max(abs(err_Pg_Qg));
    fprintf('The maximum discrepency for the Vmag and Vang results are %d and %d \n',MaxErrorVmagVang);
    fprintf('The maximum discrepency for the Pg and Qg results are %d and %d \n',MaxErrorPgQg);
end
