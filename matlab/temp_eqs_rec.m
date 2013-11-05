function [temp_dot,df_dy_temp] = temp_eqs_rec(Vr,Vi,temps,ps,opt)
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
y_tt = sparse(status.*( y_series + 1i*B/2));
y_ff = sparse(status.*( y_tt ./ tap.^2));
y_ft = sparse(status.*(-y_series ./ conj(tap_shift)));

G_ff = sparse(real(y_ff));
B_ff = sparse(imag(y_ff));
G_ft = sparse(real(y_ft));
B_ft = sparse(imag(y_ft));

% calculate temp_dot
V = Vr + 1i*Vi;
% If = Yff.*V(F) + Yft.*V(T) = (G_ff.*Vr_f - B_ff.*Vi_f + G_ft.*Vr_t - B_ft.*Vi_t)...
% + 1i * (G_ff.*Vi_f + B_ff.*Vr_f + G_ft.*Vi_t + B_ft.*Vr_t) 
If          = y_ff.*V(F) + y_ft.*V(T);
Imag_f      = abs(If);
temp_dot    = Imag_f.^2 - opt.sim.temp.K*temps; 

% calculate dImag_dVr_f, dImag_dVi_f, dImag_dVr_t, and dImag_dVi_t
If_R        = real(If);
If_I        = imag(If);
dImag_dVr_f = (If_R.*G_ff + If_I.*B_ff)./Imag_f;
dImag_dVi_f = (If_R.*(-B_ff) + If_I.*G_ff)./Imag_f;
dImag_dVr_t = (If_R.*G_ft + If_I.*B_ft)./Imag_f;
dImag_dVi_t = (If_R.*(-B_ft) + If_I.*G_ft)./Imag_f;

% calculate dFtemp_dot_dVr and dFtemp_dot_dVi
% 1st-2nd columns: derivatives for Vr_f and Vi_f
dFtemp_dot_dVr_f  = 2.*Imag_f .* dImag_dVr_f;
dFtemp_dot_dVi_f  = 2.*Imag_f .* dImag_dVi_f;

% 3rd-4th columns: derivatives for Vi_f and Vi_t
dFtemp_dot_dVr_t  = 2.*Imag_f .* dImag_dVr_t;
dFtemp_dot_dVi_t  = 2.*Imag_f .* dImag_dVi_t;


df_dy_temp          = [dFtemp_dot_dVr_f ...
                            dFtemp_dot_dVi_f ...
                                dFtemp_dot_dVr_t ...
                                    dFtemp_dot_dVi_t];

return
