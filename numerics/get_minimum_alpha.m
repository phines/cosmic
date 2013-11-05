function [ min_step_size ] = get_minimum_alpha(pointer_gx, xk, x1, x2, pk)
%evaluates a cubic interpolation of the function input. x1 and x2 being the
%gusstimates of the boundary points of minima. Outputs the minimum step
%size for minimizer.
z1 = xk + (pk.*x1);
z2 = xk + (pk.*x2);

%takes in the vector formof k1, k2
[b1,j1] = pointer_gx(z1);
[b2,j2] = pointer_gx(z2);
%[bk,jk] = pointer_gx(xk);

%converts b1, b2 and gives scalars k1, k2
k1 = (j1*pk)'*b1;
k2 = (j2*pk)'*b2;

% func_f = @(alpha1) 0.5*sum((pointer_gx(xk + alpha1*pk)).^2);
% alpha1 = 0:0.01:5;
% plot_func_f = zeros(size(alpha1));
% for i = 1:length(alpha1)
%     plot_func_f(i) = func_f(alpha1(i));
% end
%
% plot_func_f = 0.5*sum(func_f.^2);
% figure(1); clf;
% plot (alpha1,plot_func_f); hold on;


y1 = 0.5*sum(b1.^2);
y2 = 0.5*sum(b2.^2);
%yk = 0.5*sum(bk.^2);



%cubic_spline - gets a interpolated cubic spline for two given points x1 and
%x2, y1 & y2 are the values of function and k1 and k2 are vaues of
%derivatives(first) of function of the function at x1 and x2 respectively
A = [(x1).^3, (x1).^2, (x1), 1;
    (x2).^3, (x2).^2, (x2), 1;
    3.*((x1).^2), 2.*x1, 1, 0;
    3.*((x2).^2), 2.*x2, 1, 0];
%C = [c3; c2; c1; c0];
B = [y1; y2; k1; k2];

C = A\B;

% if 1
%     t = 0:.01:5;
%     cubic_spline = (c3 * ((t).^3)) + (c2 * ((t).^2)) + (c1 * (t)) + (c0);
%     %figure (2);
%     plot ((t), cubic_spline,'r');
% end

%gives out the minimum value of the solution. Solution being points where
%the derivative of the spline goes to 0.
a1 = C(2);
a2 = C(3);
a3 = C(4);


%coeffecients af the polynomial after finding the derivative.
r = 3 * a3;
s = 2 * a2;
t = 1 * a1;


solution = zeros(2,1);
d = sqrt(s^2 - 4*r*t);
solution(1) = ( -s + d ) / (2*r);
solution(2) = ( -s - d ) / (2*r);

p = solution(1);
q = solution(2);
if ~isreal(p)
    %keyboard
end

assert(isreal(p),'The solutions are complex');
assert(p>=0 || q >=0,'Both solutions are negative');
if p>q
    min_step_size = p;
else
    min_step_size = q;
end

end


