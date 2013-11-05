function opt = psoptions(opt)
% usage: opt = psoptions(opt)
% some options for the power system simulator files
if nargin==0
    opt = struct;
end

%% start with numerics options
opt = numerics_options(opt);

%% power flow options
opt.pf.tolerance = 1e-9; % convergence tolerance
opt.pf.max_iters = 20; % max power flow iterations
opt.pf.CalcIslands = 1; % iteratively calculates each island in runPowerFlow
opt.pf.CascadingPowerFlow = 0;
opt.pf.flat_start = 0;
opt.pf.load_shed_rate = 0.25; % the rate at which under frequency load shedding is done in CascadingPowerFlow mode
opt.pf.linesearch = 'backtrack';
opt.pf.update = true;

%% optimal power flow options
opt.opf.generator_commitment = 0; % switch generators on/off using MIP
opt.opf.branch_switching = 0;     % switch branches on/off using MIP

%% other options
opt.verbose = 1;
opt.seecascade = 1;

%% time-domain simulation options
opt.sim.integration_scheme = 1; % 1 = trapezoidal rule, 2 = implicit ode 15i, 2 = explicit ode15s
opt.sim.var_step = true;        % fixed (false) or variable (true) integration step size
opt.sim.ramp_frac = 0.05;       % fraction of generator allowed to ramp between generations
opt.sim.writelog  = true;       % write differential and algebraic variables to a file
opt.sim.dt_default = 1/30;      % default sampling rate (30Hz)
opt.sim.max_iters = 20;         % default number of newton iterations for trapezoidal solver
opt.sim.tolerance = 1e-6;       % newton convergence tolerance
opt.sim.draw = true;
opt.sim.overload_time_limit = 10*60; % number of seconds that the branch can sit at its rateC level (PSS/E manual)
opt.sim.t_eps = 1e-16;
% temperature relay settings
opt.sim.temp.period                 = 10;           % default is 10 minutes to reach rateA temperature
opt.sim.temp.rateA_rateB_factor     = 1.1;          % ratio between rateA and rateB limits
opt.sim.temp.case_calibrate         = 0.2127775;    % calculated previously for case2383...probably should revisit.
opt.sim.temp.K                      = ...
    -log(1-(1/opt.sim.temp.rateA_rateB_factor)^2)/(60*opt.sim.temp.period);     % branch temperature constant
% other relays' settings
opt.sim.uvls_limit = 0.5;  % threshold at which under voltage load shedding occurs
opt.sim.uvls_delta = 0.25; % the fraction of load that is shed during simulation if under voltage
opt.sim.ufls_limit = 0.01; % threshold at which under frequency load shedding occurs
opt.sim.ufls_delta = 0.25; % load shedding fraction
opt.sim.zone1_distance = 6.25; % zone 1 default distance
%opt.sim.zone1_distance = 0.9; % zone 1 default distance

% legacy
opt.simdc = opt.sim;
