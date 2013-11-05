function [ min_step_size ] = minstep(pointer_gx, xk, alphaStart, alphaEnd, pk)

z1 = xk + (pk.*alphaStart);
z2 = xk + (pk.*alphaEnd);

%takes in the vector formof k1, k2
[b1,j1] = pointer_gx(z1);
[b2,j2] = pointer_gx(z2);
%[bk,jk] = pointer_gx(xk);

%converts b1, b2 and gives scalars k1, k2
k1 = (j1*pk)'*b1;
k2 = (j2*pk)'*b2;

if ~isreal(k1) || ~isreal(k2)
    keyboard
end

func_f = @(alpha1) 0.5*sum((pointer_gx(xk + alpha1*pk)).^2);
alpha1 = 0:0.01:5;
plot_func_f = zeros(size(alpha1));
for i = 1:length(alpha1)
    plot_func_f(i) = func_f(alpha1(i));
end


%plot_func_f = 0.5*sum(func_f.^2);
figure(1); clf;
plot (alpha1,plot_func_f); hold on;


y1 = 0.5*sum(b1.^2);
y2 = 0.5*sum(b2.^2);
%yk = 0.5*sum(bk.^2);

[~, C] = evalspline(y1, y2 ,k1 ,k2, alphaStart, alphaEnd);

min_step_size = evalminspline(C);





end

