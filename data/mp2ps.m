function ps = mp2ps( mpfile )
% usage: ps = mp2ps( mpfile )
%  converts MATPOWER data to a power system (ps) structure

C 	= psconstants;
Cmp = mpConstants;
j 	= 1i;

%% get the data matrices
if  ischar(mpfile)
    mpc = feval(mpfile);
elseif isstruct(mpfile_baseMVA)
    mpc = mpfile;
else
    error('input argument does not seem to be a MATPOWER file');
end
baseMVA = mpc.baseMVA;
bus = mpc.bus;
gen = mpc.gen;
branch = mpc.branch;
if isfield(mpc,'areas'), areas = mpc.areas; end;
if isfield(mpc,'gencost'), gencost = mpc.gencost; end;

%% record the data
nBus    = size(bus,1);
nBranch = size(branch,1);
nGen    = size(gen,1);

%% record the baseMVA
ps.baseMVA = baseMVA;

%% bus
ps.bus = bus; % bus is the same
ps.bus(:,C.bu.Pd:C.bu.Bs) = 0; % zero out these values to minimize confusion
ps.bus_i = sparse( max(ps.bus(:,1)), 1 );
ps.bus_i(ps.bus(:,1)) = (1:nBus)';
Sd_bus = bus(:,Cmp.PD) + j*bus(:,Cmp.QD);
Yd_bus = bus(:,Cmp.GS) + j*bus(:,Cmp.BS);

%% branch
ps.branch = zeros(nBranch,C.br.cols);
% cols 1:11 are the same
ps.branch(:,1:11) = branch(:,1:11);
% other columns are not needed here

%% gen
ps.gen = zeros(nGen,C.ge.cols);
gcols = size(gen,2);
ps.gen(:,1:gcols) = gen(:,1:gcols);
% get gen types
ps.gen(:,C.ge.type) = bus(ps.bus_i(ps.gen(:,1)),Cmp.BUS_TYPE);
% get voltage set points
ps.gen(:,C.ge.Vsp)  = gen(:,Cmp.VG);
% status
ps.gen(:,C.ge.status) = gen(:,Cmp.GEN_STATUS);
% get the generator linear cost term
if all(gencost(:,1)==2)
    % just extract the linear cost function
    % 2-term cost functions
    subset = (gencost(:,4)==2);
    if any(subset)
        ps.gen(subset,C.ge.cost) = gencost(subset,5);
    end
    % 3-term cost functions
    subset = (gencost(:,4)==3);
    if any(subset)
        ps.gen(subset,C.ge.cost) = gencost(subset,6);
    end
    % 4-term cost functions
    subset = (gencost(:,4)==4);
    if any(subset)
        ps.gen(subset,C.ge.cost) = gencost(subset,7);
    end
else
    disp('ignoring cost data');
end

%% save the gencost and area data as is
ps.gencost = gencost;
ps.areas = areas;

%% shunt
[i_S,~,Sd] = find(Sd_bus);
[i_Y,~,Yd] = find(Yd_bus);
nShunt = length(i_S) + length(i_Y);
ps.shunt = zeros(nShunt,C.sh.cols);
ps.shunt(:,C.sh.bus) = ps.bus([i_S;i_Y],1);
ps.shunt(:,C.sh.P)   = [real(Sd);real(Yd)];
ps.shunt(:,C.sh.Q)   = [imag(Sd);imag(Yd)];
ps.shunt(:,C.sh.frac_S) = (1:nShunt)' <= length(i_S);
ps.shunt(:,C.sh.frac_Y) = not( ps.shunt(:,C.sh.frac_S) );
ps.shunt(:,C.sh.type)   = C.NO_CONTROL;
ps.shunt(:,C.sh.status) = 1;
ps.shunt(:,C.sh.value)  = C.DEFAULT_VALUE;

%% update the system
ps = updateps(ps);
