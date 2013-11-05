function [t_steps,X,Y,Z] = solve_dae(f,g,h,x0,y0,t_span,opt)
% usage: [t_steps,X,Y,Z] = solve_dae(f,g,h,x0,y0,t_span,opt)
% this function uses the first order trapezoidal rule to solve
% the dae defined by the following:
%
%  f - the differential variables
%  g - the algebraic variables
%  h - the integer variables
%  x0 - the initial set of differential variables
%  y0 - the initial set of algebraic variables
%  t_span - time span
%  opt - options

% process default inputs and outputs
if nargin<8, opt = numerics_options; end;
if isempty(h), check_discrete=false; else check_discrete=true; end
Z = [];
relay_event_hit = false;
relay_event_crossed = false;

% constants, sizes and defaults
max_newton_its  = opt.sim.max_iters; % maximum number of newton iterations
abs_tol         = opt.sim.tolerance; % newton convergence tolerance
dif_tol_min     = 0.01;              % increase step size
dif_tol_max     = 0.05;              % decrease step size
dt_min          = 0.005;             % default min step size
dt_max          = 1;                 % default max step size
alpha_0         = 1;
min_alpha       = 2^-10;
nx              = length(x0);
X               = x0;
Y               = y0;
xy0 = [x0;y0];

% find starting point
t0 = t_span(1);

% default time step
t_final = t_span(end);
dt0     = dt_max*0.10;
t_steps = t0;
% get the initial state of the dae system
f0 = f(t0,x0,y0);
if check_discrete
    z0 = h(t0,xy0);
end

% do the LU on the jacobian to get the symbolic work done
% [L,U,p,q,R] = lu(jacobian);

% set initial step size
dt = dt0;
while t0<t_final
    % choose a guess for f,g at next point
    t1 = t0 + dt;
    if t1 > t_final
        t1 = t_final;
        dt = t1-t0;
    end  % integrate up to t_final
    x1 = x0 + f0*dt;    % could use second order estimate here
    y1 = y0;            % could solve for something more intelligent
    % print something
    fprintf(' t = %f sec. and dt = %f sec...\n',t0,dt);
    
    % iteratively update x1 and y1 until it converges
    newton_it = 0;
    while newton_it < max_newton_its+1      
        % get the next residuals, jacobian blocks
        [f1,df_dx1,df_dy1] = f(t1,x1,y1);
        [g1,dg_dx1,dg_dy1] = g(t1,x1,y1);
        % check the mismatch and quit
        trapz_mismatch  = [x0 - x1 + (dt/2).*f0 + (dt/2).*f1;   % f portion
            g1];                                                % g portion
        max_trap_mis = max(abs(trapz_mismatch));
        fprintf('   Mismatch = %g\n',max_trap_mis);
        if max_trap_mis < abs_tol
            %fprintf('.converged in %d iterations with alpha=%g ...\n',newton_it,alpha);
            % if it converged, first check for any events at t0<th<t1
            % [outputs] = get_event_thresholds(inputs)
            break;
        end
        % bail out if we are at max number of iterations
        if newton_it == max_newton_its
            error(' time step did not converge...');
        end
        
        % build the jacobian that will be used in Newton's method
        trapz_df_dx = -speye(nx) + (dt/2).*df_dx1;
        trapz_df_dy = (dt/2).*df_dy1;
        jacobian    = [trapz_df_dx trapz_df_dy;
            dg_dx1 dg_dy1];
        
        % find the search direction
        p = -(jacobian\trapz_mismatch);
        % implement the backstepping to get the Newton step size
        alpha = alpha_0;
        while 1
            x1 = x1 + alpha.*p(1:nx);
            y1 = y1 + alpha.*p((nx+1):end);
            f1 = f(t1,x1,y1);
            f_mis = x0 - x1 + (dt/2).*f0 + (dt/2).*f1;
            g1 = g(t1,x1,y1);
            old_trap_mis = max_trap_mis;
            max_trap_mis = max(max(abs(f_mis)),max(abs(g1)));
            
            if max_trap_mis < old_trap_mis
                %fprintf('Newton step completed.\n');
                break
            else
                fprintf('Reducing Newton step size.\n');
                alpha = alpha/2;
            end
            if alpha<min_alpha
                fprintf('Algorithm failure in the newton step.\n');
                return
            end
        end
        newton_it = newton_it + 1;
    end
    
    % at this point assume that the solution is a good one.
    solution_good = true;
    % adjust step size if needed
    if opt.sim.var_step
        if (max(abs(f1-f0)) >= dif_tol_max) && (dt > dt_min)
            % reduce step size and go back (ie trash x1, f1, etc)
            dt = max( dt/2, dt_min);
            solution_good = false;
        else
            if max(abs(f1-f0)) < dif_tol_min
                % increase next step size
                dt = min( 2*dt, dt_max);
            end
        end
    end
    % check to see if a threshold was crossed
    if check_discrete
        xy1 = [x1;y1];
        z1 = h(t1,xy1);
        % check to see if we hit a threshold
        hit_thresh = abs(z1)<opt.sim.eps_thresh;
        crossed = (z1<= -opt.sim.eps_thresh);
        if any(hit_thresh)
            relay_event_hit = true;
        end
        % figure out when the threshold crossed
        if any(crossed) && dt > 1/100
            dz = z1(crossed) - z0(crossed);
            t_thresh = -(dt./dz).*z0(crossed) + t0;
            [~,which_thresh] = min(t_thresh);
            dt_temp = t_thresh(which_thresh) - t0;
            dt = min(dt,dt_temp);
            solution_good = false;
        elseif any(crossed) && dt <= 1/100
            z0 = z1;
            relay_event_crossed = true;
        end       
    end
    % record the solution for this step
    if solution_good
        if relay_event_hit || relay_event_crossed
            relay_event = or(hit_thresh,crossed);
        else
            relay_event = [];
        end
        % commit the results to memory
        % if it converged and no events, then save t1, x1 and y1; advance
        t0          = t1;           % commit time advance
        x0          = x1;
        y0          = y1;
        f0          = f1;
        % save the variables
        X           = [X x1];       %#ok<*AGROW>
        Y           = [Y y1];
        t_steps     = [t_steps;t1];
        if check_discrete && any(relay_event)
            keyboard
            Z = z1;
            break;
        end
    end
end
