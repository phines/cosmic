function [g,dg_dx] = mismatch_rec(x,Ybus,Vmag,Vr,Vi,Sg,Sd,pq,pv,ref,PartFact)
% usage: [g,dg_dx] = rec_mismatch(x,Ybus,Vmag,Sg,Sd,pq,pv,ref,PartFact)
%
% a power flow mismatch function, with participation factors
% the buses that are not in pq and pv are assumed to be generators
% that do load-following. The first of these will get the
% zero angle reference. The bus at "ref" will be assiged a zero
% angle
% If PartFact (participation factors) is defined it should be
% an n x 1 vector indicating the extent to which the generators particpate
% in load following

% constants
j = 1i;
nBus = size(Ybus,1);

if nargin<9
    PartFact = ones(nBus,1);
end
% convert pv/pq to logical arrays
if length(pv)<nBus
    pv_ = false(nBus,1);
    pv_(pv) = true;
    pv = pv_;
end
if length(pq)<nBus
    pq_ = false(nBus,1);
    pq_(pq) = true;
    pq = pq_;
end
if length(ref)<nBus
    ref_ = false(nBus,1);
    ref_(ref) = true;
    ref = ref_;
end
if (sum(ref)~=1), error('Must have only one ref bus'); end

% build an index
ix.Vr = 1:(nBus-1);
ix.Vi = (nBus-1) + (1:(nBus-1));
ix.rho = max(ix.Vi) + 1; 

% degenerate case: 
if isempty(ix.rho) 
    ix.rho = max(ix.Vr) + 1;
end

nx = 2.*nBus - 1 ; 

% extract things from x
Vr = sparse(Vr);
Vi = sparse(Vi);
Vr(~ref) = x(ix.Vr);
Vi(~ref) = x(ix.Vi);
rho = x(ix.rho);
Vmag(pq) = sqrt(Vr(pq).^2+Vi(pq).^2);


% calculate the voltage
V = Vr + j*Vi; 
% make adjustments to Sg for gen ramping
ramp_gen        = ~pv & ~pq;
% calculate the total load according to ZIPE matrix
zipe_cols       = size(Sd,2);
if zipe_cols == 1
    S_zipe = Sd;
elseif zipe_cols == 5
    S_Z = Sd(:,1) .* Vmag.^2;
    S_I = Sd(:,2) .* Vmag;
    S_P = Sd(:,3);
    S_E = Sd(:,4) .* Vmag.^Sd(:,5);
    S_zipe = S_Z + S_I + S_P + S_E;
else
    error('zipe load model matrix is not the right size');
end

Sbus            = Sg - S_zipe;
Sbus(ramp_gen)  = Sbus(ramp_gen) + PartFact(ramp_gen)*rho;
% compute the final mismatch
miscx = (V .* conj(Ybus * V)) - (Sbus);  

% g = [real(miscx); imag(miscx(pq))];
g = [real(miscx); imag(miscx(~ref))];

% find the locations for Q(pv) in mismatch vector, then replace them with 
% (delV)^2 = (Vr(pv).^2 + Vi(pv).^2) - Vmag(pv).^2;
pv_loc = [false(nBus,1); pv(~ref)];
g(pv_loc) = (Vr(pv).^2 + Vi(pv).^2) - Vmag(pv).^2;

% Jacobian
if nargout>1
    n = nBus; 
    
    % do some matrix algebra 
    G = real(Ybus);
    B = imag(Ybus);
    diagVr     = spdiags(Vr, 0, n, n);
    diagVi     = spdiags(Vi, 0, n, n);
    
    % computes dP/dVr, dp/dVi, dQ/Vr, and dQ/Vi
    dP_dVr = diagVr*G + diag(G*Vr) - diag(B*Vi) + diagVi*B;
    dP_dVi = -diagVr*B + diag(B*Vr) + diag(G*Vi) + diagVi*G;
    dQ_dVr = diagVi*G - diag(G*Vi) - diag(B*Vr) - diagVr*B;
    dQ_dVi = -diagVr*G + diag(G*Vr) - diag(B*Vi) - diagVi*B;

    % now put these into dg_dx
    dg_dx = sparse(nx,nx);
    dP_rows = (1:n);
    %dQ_rows = (1:npq) + n;
    
    % dP_dVr
    [rows,cols,values] = find(dP_dVr(:,~ref));
    dg_dx = dg_dx + sparse(rows,cols,values,nx,nx);

    % dP_dVi
    [rows,cols,values] = find(dP_dVi(:,~ref));
    dg_dx = dg_dx + sparse(rows,cols+(nBus-1),values,nx,nx); 

    % dQ_dVr
    dQ_dVr(pv,:) = 0;
    PVrow = find(pv);
    PVcol = PVrow';
    PVvalue_Vr = 2*Vr(pv);
    dQ_dVr = dQ_dVr + sparse(PVrow,PVcol,PVvalue_Vr,n,n);
    [rows,cols,values] = find(dQ_dVr(~ref,~ref));
    dg_dx = dg_dx + sparse(rows+nBus,cols,values,nx,nx);

    % dQ_dVi
    dQ_dVi(pv,:) = 0;
    PVvalue_Vi = 2*Vi(pv);
    dQ_dVi = dQ_dVi + sparse(PVrow,PVcol,PVvalue_Vi,n,n);    
    [rows,cols,values] = find(dQ_dVi(~ref,~ref));
    dg_dx = dg_dx + sparse(rows+nBus,cols+(nBus-1),values,nx,nx);  
    
    % dP_drho
    dg_dx = dg_dx + sparse(dP_rows(ramp_gen),ix.rho,-PartFact(ramp_gen),nx,nx);
   
    % fix the derivatives with ZIP[E] contributions
    dP_E_dVr = real(Sd(:,4)) .* Sd(:,5) .* Vr .* Vmag.^(2.*(Sd(:,5)./2-1));
    dP_E_dVi = real(Sd(:,4)) .* Sd(:,5) .* Vi .* Vmag.^(2.*(Sd(:,5)./2-1));
    dQ_E_dVr = imag(Sd(:,4)) .* Sd(:,5) .* Vr .* Vmag.^(2.*(Sd(:,5)./2-1));
    dQ_E_dVi = imag(Sd(:,4)) .* Sd(:,5) .* Vi .* Vmag.^(2.*(Sd(:,5)./2-1));
    % fix the derivatives with [Z]IPE contributions
    dP_Z_dVr = 2 .* real(Sd(:,1)) .* Vr;
    dP_Z_dVi = 2 .* real(Sd(:,1)) .* Vi;
    dQ_Z_dVr = 2 .* imag(Sd(:,1)) .* Vr;
    dQ_Z_dVi = 2 .* imag(Sd(:,1)) .* Vi;
    % fix the derivatives with Z[I]PE contributions
    dP_I_dVr = real(Sd(:,2)) .* Vr ./Vmag;
    dP_I_dVi = real(Sd(:,2)) .* Vi ./Vmag;
    dQ_I_dVr = imag(Sd(:,2)) .* Vr ./Vmag;
    dQ_I_dVi = imag(Sd(:,2)) .* Vi ./Vmag;
    
    SS_zipe = S_zipe; 
    SS_zipe(~pq)= 0;

    rows = find(SS_zipe); cols = find(SS_zipe(~ref));
    dg_dx = dg_dx + sparse(rows,cols,dP_E_dVr(rows),nx,nx);
    dg_dx = dg_dx + sparse(rows,cols+(nBus-1),dP_E_dVi(rows),nx,nx);
    dg_dx = dg_dx + sparse(rows,cols,dP_Z_dVr(rows),nx,nx);
    dg_dx = dg_dx + sparse(rows,cols+(nBus-1),dP_Z_dVi(rows),nx,nx);
    dg_dx = dg_dx + sparse(rows,cols,dP_I_dVr(rows),nx,nx);
    dg_dx = dg_dx + sparse(rows,cols+(nBus-1),dP_I_dVi(rows),nx,nx);
    
    rows = find(SS_zipe(~ref)); cols = rows;
    dg_dx = dg_dx + sparse(rows+nBus,cols,dQ_E_dVr(rows),nx,nx);
    dg_dx = dg_dx + sparse(rows+nBus,cols+(nBus-1),dQ_E_dVi(rows),nx,nx);    
    dg_dx = dg_dx + sparse(rows+nBus,cols,dQ_Z_dVr(rows),nx,nx);
    dg_dx = dg_dx + sparse(rows+nBus,cols+(nBus-1),dQ_Z_dVi(rows),nx,nx);  
    dg_dx = dg_dx + sparse(rows+nBus,cols,dQ_I_dVr(rows),nx,nx);
    dg_dx = dg_dx + sparse(rows+nBus,cols+(nBus-1),dQ_I_dVi(rows),nx,nx);  
 
end

