function [ps,discrete] = process_event(ps,event,opt)
% usage: [ps,discrete] = process_event(ps,event,opt)
% process the event specified in the event vector (one row of an event matrix)


if size(event,1)>1 
    error('process_event:err','can only process one event at a time.');
end

verbose = opt.verbose;
discrete = false;   % is this a discrete event?

% extract some data
C = psconstants;
t = event(C.ev.time);

% record the event
ps.event_record = [ps.event_record;event];

% do the processing
switch event(C.ev.type)
    case C.ev.start
        if verbose, fprintf(' simulation start event.'); end

    case C.ev.fault     % three phase to ground fault by default
        error('process_event:err','faults not supported yet.');
        
    case C.ev.trip_branch
        branch_id = event(C.ev.branch_loc);
        branch_ix = ps.branch_i(branch_id);
        if ps.branch(branch_ix,C.br.status) ~= 0
            ps.branch(branch_ix,C.br.status) = 0;
            discrete = true;
            % output the branch id
            if verbose, fprintf('  t = %.2f: Branch %d tripped...\n',t,branch_id); end
        end
        
    case C.ev.close_branch
        branch_id = event(C.ev.branch_loc);
        branch_ix = ps.branch_i(branch_id);
        if ps.branch(branch_ix,C.br.status) ~= 1
            ps.branch(branch_ix,C.br.status) = 1;
            discrete = true;
            % output the branch id
            if verbose, fprintf('  t = %.2f: Branch %d closed...\n',t,branch_id); end
        end
            
    case C.ev.trip_bus
        bus_no = event(1,C.ev.bus_loc);
        br_set = any( ps.branch(:,1:2)==bus_no, 2 );
        ps.branch(br_set,C.br.status) = 0;
        % trip gens and shunts at this bus
        ps.gen(ps.gen(:,1)==bus_no,C.gen.status) = 0;
        ps.shunt(ps.shunt(:,1)==bus_no,C.shunt.status) = 0;
        discrete = true;
        if verbose, fprintf('  t = %.2f: Bus %d tripped...\n',t,bus_no); end
        
    case C.ev.trip_gen
        gen_id = event(C.ev.gen_loc);
        gen_ix = ps.gen_i(gen_id);
        if ps.gen(gen_ix,C.gen.status)~=0
            ps.gen(gen_ix,C.gen.status) = 0;
            ps.gen(gen_ix,C.Pg) = 0;
            ps.gen(gen_ix,C.Qg) = 0;
            if verbose, fprintf('  t = %.2f: Gen %d tripped...\n',t,gen_id); end
            discrete = true;
        end
        
    case C.ev.shed_load
        shunt_id = event(C.ev.shunt_loc);
        shunt_ix = ps.shunt_i(shunt_id);
        prev_load = ps.shunt(shunt_ix,C.sh.P).*ps.shunt(shunt_ix,C.sh.factor); 
        ps.shunt(shunt_ix,C.sh.factor) = max(ps.shunt(shunt_ix,C.sh.factor) - event(C.ev.quantity), 0);
        curr_load = ps.shunt(shunt_ix,C.sh.P).*ps.shunt(shunt_ix,C.sh.factor); 
        discrete = true;
        bus_no = ps.shunt(shunt_ix,1);
        if opt.verbose, fprintf('  t = %.2f: %.2f MW of load shedding at bus %d...\n',t,prev_load-curr_load,bus_no); end

    case C.ev.relay_trigger
        re_index = event(C.ev.relay_loc);
        re_type  = ps.relay(re_index,C.re.type);
        switch re_type
            case C.re.dist
                % 1) zone 1 distance relay
                branch_id = ps.relay(re_index,C.re.branch_loc);
                ps.relay(re_index,C.re.tripped) = 1;
                % create the branch trip event
                new_event = zeros(1,C.ev.cols);
                new_event(C.ev.time) = t;
                new_event(C.ev.type) = C.ev.trip_branch;
                new_event(C.ev.branch_loc) = branch_id;
                discrete = true;
                % tell the user
                if verbose, fprintf('  t = %.2f: Distance relay trip at branch %d...\n',t,branch_id); end
                % process the event
                ps = process_event(ps,new_event,opt);
            case C.re.ufls
                % 2) under frequency load shedding
                % TODO--figure out how to implement this
            case C.re.uvls
                % 3) undervoltage load shedding
                shunt_id = ps.relay(re_index,C.re.shunt_loc);
                shunt_ix = ps.shunt_i(shunt_id);
                bus_id   = ps.shunt(ps.shunt_i(shunt_id),1);
                % create the load shedding event
                new_event = zeros(1,C.ev.cols);
                new_event(C.ev.time) = t;
                new_event(C.ev.type) = C.ev.shed_load;
                new_event(C.ev.shunt_loc) = shunt_id;
                new_event(C.ev.quantity)  = ps.relay(re_index,C.re.setting1);
                discrete = true;
                % flag the relay as tripped
                if ps.shunt(shunt_ix,C.sh.factor)<=new_event(C.ev.quantity)
                    ps.relay(re_index,C.re.tripped) = 1;
                end
                % tell user:
                if opt.verbose, fprintf('  t = %.2f: Undervoltage load shedding triggered at bus %d...\n',t,bus_id); end
                % process the event
                ps = process_event(ps,new_event,opt);
            case C.re.temp
                % 4) temp relay
                branch_id = ps.relay(re_index,C.re.branch_loc);
                % create the branch trip event
                new_event = zeros(1,C.ev.cols);
                new_event(C.ev.time) = t;
                new_event(C.ev.type) = C.ev.trip_branch;
                new_event(C.ev.branch_loc) = branch_id;
                discrete = true;
                ps.relay(re_index,C.re.tripped) = 1;
                % tell the user
                if verbose, fprintf('  t = %.2f: Temperature relay trip at branch %d...\n',t,branch_id); end
                % process the event
                ps = process_event(ps,new_event,opt);
        end
    otherwise
        error('process_event:err','Unknown event type');
end
