function fg = differential_algebraic_rec(t,xy,nx,ps,opt)
% usage: fg = differential_algebraic_rec(t,xy,nx,ps,opt)

x = xy(1:nx);
y = xy((nx+1):end);

fg = [differential_eqs_rec(t,x,y,ps,opt);
         algebraic_eqs_rec(t,x,y,ps,opt);];
