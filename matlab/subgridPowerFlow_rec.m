function [Vmag,theta,Vr,Vi,success] = subgridPowerFlow_rec(bus_subset,Ybus,Vmag,theta,Vr,Vi,Sg,Sd,pq,pv,ref,part_fact,nr_opts)
% usage: [Vmag,theta,success] = subgridPowerFlow(bus_subset,Ybus,Vmag,theta,Sg,Sd,pq,pv,ref,part_fact,nr_opts)
% power flow for a subgrid defined by bus_subset
n_sub = sum(bus_subset);

if sum(ref)==0
    Vmag(bus_subset)  = 0;
    theta(bus_subset) = 0;
    Vr(bus_subset)  = 0;
    Vi(bus_subset)  = 0;
    success = false;
elseif n_sub==1                  % we have one generator bus
    theta(bus_subset) = 0;
    Vi(bus_subset) = 0;
    success = true;
else
    % build the decision vector
    ix = struct;
    nBus = size(Ybus,1);
    ix.Vr = 1:(nBus-1);
    ix.Vi = (nBus-1) + (1:(nBus-1));
    ix.rho = max(ix.Vi) + 1;     % generator ramping variable

    if isempty(ix.rho)           % degenerate case
        ix.rho = 1 + max(ix.Vr);
    end    
    nx = max(ix.rho);
    x = zeros(nx,1);
    x(ix.Vr)   = Vr(~ref & bus_subset);
    x(ix.Vi)   = Vi(~ref & bus_subset);
    x(ix.rho)  = 0;

    if any(x(ix.Vr))==0          % missing initial guess
        x(ix.Vr) = 1;
        x(ix.Vi) = 0;
    end

    % subset the variables for this subgrid
    Ybus_sub    = Ybus(bus_subset,bus_subset);
    Vmag_sub    = Vmag(bus_subset);
    Vr_sub      = Vr(bus_subset);
    Vi_sub      = Vi(bus_subset);
    Sg_sub      = Sg(bus_subset);
    Sd_sub      = Sd(bus_subset,:);
    pq_sub      = pq(bus_subset);
    pv_sub      = pv(bus_subset);
    ref_sub     = ref(bus_subset);
    
    if part_fact(ref)==0
        part_fact(ref) = 1;
    end

    % set up a virtual function to solve
    g = @(newx)mismatch_rec(newx,Ybus_sub,Vmag_sub,Vr_sub,Vi_sub,Sg_sub,Sd_sub,pq_sub,pv_sub,ref_sub,part_fact);
    % try to solve the power flow problem
    [x,success] = nrsolve(g,x,nr_opts);
    
    % if this failed re-try with a flat start
    if ~success
        if nr_opts.nr.verbose
           disp('Re-trying power-flow with a flat start');
        end
        x(ix.Vi) = 0;
        x(ix.Vr) = 1;
        [x,success] = nrsolve(g,x,nr_opts);
    end
    % save the results
    Vr(~ref & bus_subset) = x(ix.Vr);
    Vi(~ref & bus_subset) = x(ix.Vi);
    Vi(ref) = 0;
end