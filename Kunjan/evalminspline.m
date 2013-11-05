function [ min_value ] = evalminspline(C)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%c0 = C(1);
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
    keyboard
end

assert(isreal(p),'The solutions are complex');
assert(p>=0 || q >=0,'Both solutions are negative');
if p>q
    min_value = p;
else
    min_value = q;
end

return

