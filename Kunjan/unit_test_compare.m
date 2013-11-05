function unit_test_compare(ps,fname1,fname2,casename)
% Compares the cosmic simulation results in two files to see if they look similar

% get from fname1
[t,delta,omega,~,~,~,Vmag,theta] = read_outfile(fname1,ps);

% get from fname2
[t_u,delta_u,omega_u,~,~,~,Vmag_u,theta_u] = read_outfile(fname2,ps);

% put the unit test file results onto the calculated time scale
delta_u = interp1(t_u,delta_u,t,'pchip');
omega_u = interp1(t_u,omega_u,t,'pchip');
Vmag_u  = interp1(t_u,Vmag_u ,t,'pchip');
theta_u = interp1(t_u,theta_u,t,'pchip');

%% check the results
success=true;
if ~all((delta(:) - delta_u(:)) < 1e-6)
    fprintf('Error in %s: found discrepancy on machine angles \n',casename);
    success=false;
end
if ~all((omega(:) - omega_u(:)) < 1e-6)
    fprintf('Error in %s: found discrepancy on machine speeds \n',casename);
    success=false;
end
if ~all((Vmag(:) - Vmag_u(:)) < 1e-6)
    fprintf('Error in %s: found discrepancy on voltage magnitudes \n',casename);
    success=false;
end
if ~all((theta(:) - theta_u(:)) < 1e-6)
    fprintf('Error in %s: found discrepancy on voltage angles \n',casename);
    success=false;
end
if success
    fprintf('----------------------------------------------------\n\n');
    fprintf('Successfully verified %s... \n\n',casename);
else
    fprintf('Unit test failed\n');
    keyboard
end

fprintf('Unit test for %s complete\n',casename);
fprintf('----------------------------------------------------\n\n');
