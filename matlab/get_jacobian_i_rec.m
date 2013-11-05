function [dF_dxy,dF_dxyp] = get_jacobian_i_rec(t,xy,xyp,ps,opt)
% usage: [df_dxy,dF_dxyp] = get_jacobian_i_rec(t,xy,xyp,ps,opt)

[~, dF_dxy,dF_dxyp] = differential_algebraic_i_rec(t,xy,xyp,ps,opt);
