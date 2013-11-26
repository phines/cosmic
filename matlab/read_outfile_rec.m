function [t,delta,omega,Pm,Eap,Vmag,theta,Vr,Vi,E1,Efd] = read_outfile_rec(fname,ps,opt)
% usage: [t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd] = read_outfile(fname,ps)

n  = size(ps.bus,1);
ng = size(ps.gen,1);
m  = size(ps.branch,1);
n_sh = size(ps.shunt,1);
ix   = get_indices_rec(n,ng,m,n_sh,opt);
j = 1i;


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

% Y vars
Vr = Y(:,ix.y.Vr);
Vi = Y(:,ix.y.Vi);
Vmag  = abs(Vr + j.*Vi);
theta = angle(Vr + j.*Vi);

