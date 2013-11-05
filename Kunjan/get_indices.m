function index = get_indices(n_bus,n_macs,n_branches,n_shunts,interleave)
% usage: index = get_indices(n_bus,n_macs,n_branches,n_shunts,interleave)
% produces a structure that helps us to know where stuff is in x,y,f,g

if nargin < 5, interleave = true; end     % debugging flag to test different reorderings of the x and y variables

% differential (generator,exciter,governor,relays) variable index
nxsmac = 6;				% number of differential variables per machine 
if interleave
    index.x.delta       = (1:nxsmac:nxsmac*n_macs);
    index.x.omega_pu    = (2:nxsmac:nxsmac*n_macs);
    index.x.Pm          = (3:nxsmac:nxsmac*n_macs);
    index.x.Eap         = (4:nxsmac:nxsmac*n_macs);
	index.x.E1 			= (5:nxsmac:nxsmac*n_macs);
	index.x.Efd			= (6:nxsmac:nxsmac*n_macs);    
    index.x.temp        = (1:n_branches) + nxsmac*n_macs;
else
    index.x.delta       = (1:n_macs);
    index.x.omega_pu    = (1:n_macs) + n_macs;
    index.x.Pm          = (1:n_macs) + n_macs*2;
    index.x.Eap         = (1:n_macs) + n_macs*3; 
	index.x.E1			= (1:n_macs) + n_macs*4;
	index.x.Efd 		= (1:n_macs) + n_macs*5; 
    index.x.temp        = (1:n_branches) + n_macs*6;
end

index.nx            = nxsmac*n_macs + n_branches;
index.x.omega       = index.x.omega_pu; 

% differential equation index is the same as x index
index.f.delta_dot = index.x.delta;
index.f.omega_dot = index.x.omega_pu;
index.f.Pm_dot    = index.x.Pm;
index.f.Eap_dot   = index.x.Eap;
index.f.E1_dot 	  = index.x.E1;
index.f.Efd_dot   = index.x.Efd;
index.f.temp_dot  = index.x.temp;

index.nf = index.nx;
index.f.swing = index.f.omega_dot;

% algebraic variable index 
if interleave
    index.y.Vmag    = (1:2:2*n_bus);
    index.y.theta   = (2:2:2*n_bus); 
else
    index.y.Vmag    = (1:n_bus);
    index.y.theta   = (1:n_bus) + n_bus;
end
index.y.delta_sys = n_bus*2 + 1;
index.ny = 2*n_bus + 1;

% algebraic function index
index.g.P       = index.y.Vmag;
index.g.Q       = index.y.theta;
index.g.slack   = index.y.delta_sys;

index.ng = index.ny;

% relay index
index.re.temp = 1:n_branches;
index.re.uvls = (1:n_shunts) + n_branches;
index.re.ufls = (1:n_macs)   + n_branches + n_shunts;
index.re.dist = (1:n_branches) + n_branches + n_shunts + n_macs;
return
