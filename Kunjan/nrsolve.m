function [x,success,k] = nrsolve(eval_g,x0,opts)
% use the Newton-Raphson method to solve a set of non-linear equations
% usage: [x,success,k] = nrsolve(g,x0,opts)

% process inputs
if nargin<2, error('at least 2 inputs needed'); end
if nargin<3, opts = numerics_options; end;
if nargout(eval_g)==1
    error('The input function needs to output the jacobian');
end

max_iters = opts.nr.max_iterations;
tolerance = opts.nr.tolerance;
verbose   = opts.nr.verbose;
use_fsolve = opts.nr.use_fsolve;
alpha = 0;

x = x0;
success = false;
tic;

%{
% check the condition of the jac
[~,J] = eval_g(x0);
c = condest(J);
if c>1e8
    warning('nrsolve:JCondition','J appears to be ill conditioned: %g',c);
end
%}

% solve with fsolve
if use_fsolve
    fsolve_opts = optimset( 'Jacobian','on',...
                            'Algorithm','trust-region-dogleg',...
                            'Display','off');
    [x,~,flag] = fsolve(eval_g,x0,fsolve_opts);
    success = (flag==1);
else
    % print something
    if verbose
        disp('Iter  Mismatch  alpha');
    end

    for k = 0:max_iters
        % evaluate the function and the derivative
        [g,J] = feval(eval_g,x);
        max_mismatch  = max(abs(g));
        mean_mismatch = mean(abs(g));
        % print something
        if verbose
            fprintf('%4d %10.7f %10.7f %g\n',k,max_mismatch,mean_mismatch,alpha);
        end
        % check for convergence
        if max_mismatch<tolerance %checks if max is is greater than tolerance
            success = true;
            if verbose
                et = toc;
                fprintf('Solution found in %g seconds\n',et);
            end
            break;
        end
        % choose the search direction
        p = -(J\g);
        % do some sort of line search to select the step size and the new x
        [x,alpha] = linesearch(p,x,eval_g,opts);
        if alpha < 1e-12
            success = false;
            break
        end
    end
end

if ~success && verbose
    disp(' Did not find a solution for the algebraic variables in a subgrid');
end

return;

function [x,alpha,g] = linesearch(p,x_prev,eval_g,opts)
% perform a line search to choose a step size

% get the starting point
f_prev = qnorm(feval(eval_g,x_prev));
% choose the method

switch opts.nr.linesearch
    case 'backtrack'
        alpha = 1;
        alpha_min = opts.nr.alpha_min;
        mu        = opts.nr.mu;
        while alpha > alpha_min
            x = x_prev + alpha*p;
            g = feval(eval_g,x);
            f = qnorm(g);
            % test the sufficient decrease (Armijo) condition
            if f <= f_prev + mu * alpha * sum( p.*x );
                break
            end
            alpha = alpha / 2;
            if 0
                a = -1:.01:1;
                for i = 1:length(a)
                    fi(i) = qnorm(feval(eval_g,x_prev + a(i)*p));
                end
                figure(1);
                plot(a,fi);
                pause
            end
        end
    case 'exact'
        f_mis = @(a) qnorm(feval(eval_g,x_prev + a*p));
        [alpha,~,flag] = fminunc(f_mis,1,optimset('Display','off','LargeScale','off'));
        if flag==0
            alpha = 0;
        end
        x = x_prev + alpha*p;
        g = feval(eval_g,x);
    case 'cubic_spline'
        [alpha] = minstep(eval_g, x_prev, 0, 1.5, p);
        x = x_prev + alpha*p;
        g = feval(eval_g,x);
        
    otherwise
        error('unsupported line search method');
end

function f = qnorm(x)
% simple quadratic norm:
% note that the derivative of this is 2*x/2 = x

f = sum( x.^2 ) / 2;

