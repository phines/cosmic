function unix_exit_flag = cosmic(case_data,event_matrix)
% This is the entry function for the cosmic simulator
% Inputs:
%  case_data gives the power system network data (see psconstants for structure information)
%  event_matrix is a matrix describing the exogenous events that trigger stuff
% Either input can be provided as a matlab matrix, or as the name of a mat file

%% read the input data
% case data
if ischar(case_data)
    load(case_data);
else
    ps = case_data;
end
ps = updateps(ps);
% event matrix
if ischar(event_matrix)
    load(event_matrix);
else
    event = event_matrix;
end

%% prepare the output
unix_exit_flag = 1;

%% prepare options
opt = psoptions;
opt.sim.dt_default = 0.1;
opt.nr.use_fsolve = true;
opt.verbose = true;
opt.writelog = false;

%% prepare the data
% run power flow
ps = newpf(ps,opt);
% get the Ybus
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
% build the machine variables
[ps.mac,ps.exc,ps.gov] 		= get_mac_state(ps,'salient');
% initialize relays
ps.relay    				= get_relays(ps);

%% run the simulation
outputs = simgrid(ps,event,'sim_case9',opt);

%% write the endogenous event file
outfilename = sprintf('events_%s_%s.m',ps.casename,outputs.start_time);
write_event_matrix(outputs.endo_event,outfilename)

%% return the exit flag
if outputs.success
    fprintf('Simulation succeeded. %.2f MW of demand lost\n',outputs.demand_lost);
    unix_exit_flag = 0;
else
    fprintf('Simulation failed.\n');
end
