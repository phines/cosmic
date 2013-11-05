function [ cubic_spline, C ] = evalspline ( y1, y2, k1, k2, x1, x2 )
%cubic_spline - gets a interpolated cubic spline for 2 given points x1 and
%x2
%   Detailed explanation goes here

%y1 = fx(x1);
%y2 = fx(x2);
%k1 = df_dx(x1);
%k2 = df_fx(x2);



A = [(x1).^3, (x1).^2, (x1), 1;
    (x2).^3, (x2).^2, (x2), 1; 
    3.*((x1).^2), 2.*x1, 1, 0;
    3.*((x2).^2), 2.*x2, 1, 0];
%C = [c3; c2; c1; c0];
B = [y1; y2; k1; k2];

C = A\B;
c0 = C(4);
c1 = C(3);
c2 = C(2);
c3 = C(1);
%cubic_spline = poly ([c0, c1, c2, c3], [t]);

if 1
    t = 0:.01:5;
    cubic_spline = (c3 * ((t).^3)) + (c2 * ((t).^2)) + (c1 * (t)) + (c0); 
    %figure (2);
    plot ((t), cubic_spline,'r');
end

return
