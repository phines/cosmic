function [temp_dot,df_dy_temp] = temp_eqs(Vmags,Thetas,temps,ps,opt)
% usage: [temp_dot,df_dy_temp] = temp_eqs(Vmags,Thetas,temps,ps,opt)
%
% This function outputs the differential temperature values and 
% derivatives wrt to the algebraic variables.
C   = psconstants;
eps = 1e-9;

% calculate dFtemp_dot_dVmag_values and dFtemp_dot_dtheta_values
n       = size(ps.bus,1); 
m       = size(ps.branch,1);
F       = ps.bus_i(ps.branch(:,C.br.from));
T       = ps.bus_i(ps.branch(:,C.br.to));
status  = (ps.branch(:,C.br.status)>=1);
R       = ps.branch(:,C.br.R);
X       = ps.branch(:,C.br.X);
B       = ps.branch(:,C.br.B);
shift   = ps.branch(:,C.br.shift);
tap     = ps.branch(:,C.br.tap);
tap( abs(tap)<eps ) = 1.0;

y_series = 1./(R+1i*X);
tap_shift = tap .* exp(-1i*pi/180 * shift);
y_tt = status.*( y_series + 1i*B/2);
y_ff = status.*( y_tt ./ tap.^2);
y_ft = status.*(-y_series ./ conj(tap_shift));

bi   = [(1:m)';(1:m)'];
Yf  = sparse(bi, [F; T], [y_ff; y_ft], m, n);

% calculate temp_dot
V           = Vmags .* exp(1i*Thetas);
If          = Yf * V;
Imag_f      = abs(If);
temp_dot    = Imag_f.^2 - opt.sim.temp.K*temps; 

%
% get info from the Ybus and form pairs of vectors [a,b]
Vmag_f      = Vmags(F);
Vmag_t      = Vmags(T);
Thetas_f    = Thetas(F);
Thetas_t    = Thetas(T);
Thetas_ft   = Thetas_f - Thetas_t;

% calculate dFtemp_dot_dVmag and dFtemp_dot_dtheta
len_a = abs(y_ff).*Vmag_f;
len_b = abs(y_ft).*Vmag_t;
len_c = Imag_f;
psi   = angle(y_ff) - angle(y_ft) + Thetas_ft;
ang_C = unwrap( pi - psi );

dImag_da    = (len_a - len_b.*cos(ang_C))./len_c;
dImag_db    = (len_b - len_a.*cos(ang_C))./len_c;
dImag_dC    = len_a.*len_b.*sin(ang_C)./len_c;

dImag_dVmag_f   = dImag_da .* abs(y_ff);
dImag_dVmag_t   = dImag_db .* abs(y_ft);

% 1st column: derivatives for Vmag_f
dFtemp_dot_dVmag_f  = 2.*Imag_f .* dImag_dVmag_f;

% 2nd column: derivatives for Vmag_t
dFtemp_dot_dVmag_t  = 2.*Imag_f .* dImag_dVmag_t;
 
% 3rd column: derivatives for theta
dFtemp_dot_dtheta   = - 2 .* Imag_f .* dImag_dC;

df_dy_temp          = [dFtemp_dot_dVmag_f ...
                            dFtemp_dot_dVmag_t ...
                                dFtemp_dot_dtheta];

return


function angle = unwrap(angle)

while angle > pi
    angle = angle - 2*pi;
end

while angle < -pi
    angle = angle + 2*pi;
end

return
%