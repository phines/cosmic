function [t,delta,omega,Pm,Eap,temp,Vmag,theta,E1,Efd] = read_outfile(fname,ps)
% usage: [t,delta,omega,Pm,Eap,temp,Vmag,theta,E1,Efd] = read_outfile(fname,ps)

n  = size(ps.bus,1);
ng = size(ps.gen,1);
m  = size(ps.branch,1);
n_sh = size(ps.shunt,1);
ix   = get_indices(n,ng,m,n_sh);


data = csvread(fname,1);

t = data(:,1);
X = data(:,1+(1:ix.nx));
Y = data(:,1+ix.nx+(1:ix.ny));

% X vars
delta = X(:,ix.x.delta);
omega = X(:,ix.x.omega)*2*pi*ps.frequency;
Pm    = X(:,ix.x.Pm);
Eap   = X(:,ix.x.Eap);
E1    = X(:,ix.x.E1);
Efd   = X(:,ix.x.Efd);
temp  = X(:,ix.x.temp);

% Y vars
Vmag  = Y(:,ix.y.Vmag);
theta = Y(:,ix.y.theta);

