function [g,dg_dx] = mismatch(x,Ybus,Vmag,Sg,Sd,pq,pv,ref,PartFact)
% usage: [g,dg_dx] = mismatch(x,Ybus,Vmag,Sg,Sd,pq,pv,ref,PartFact)
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
npq = sum(pq);
ix.theta = 1:(nBus-1);
ix.Vmag  = (1:npq) + max(ix.theta);
% add in a variable for the generator ramping
ix.rho = 1 + max(ix.Vmag);
% degenerate case:
if isempty(ix.rho) % degenerate case
    ix.rho = 1 + max(ix.theta);
end
nx = nBus + npq;

% extract things from x
theta = zeros(nBus,1);
theta(~ref) = x(ix.theta);
Vmag(pq)    = x(ix.Vmag);
rho         = x(ix.rho);

% calculate the voltage
V = Vmag.*exp(j*theta);
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
g = [real(miscx); imag(miscx(pq))];

% Jacobian
if nargout>1
    n = nBus;
    
    % do some matrix algebra (borrowed from MATPOWER)
    Ibus = Ybus * V;
    diagV     = spdiags(V, 0, n, n);
    diagIbus  = spdiags(Ibus, 0, n, n);
    diagVnorm = spdiags(exp(j*theta), 0, n, n);
    dSbus_dVm = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm;
    dSbus_dVa = j * diagV * conj(diagIbus - Ybus * diagV);
    % now put these into dg_dx
    dg_dx = sparse(nx,nx);
    dP_rows = (1:n);
    %dQ_rows = (1:npq) + n;
    % dP_dtheta
    [rows,cols,values] = find(real(dSbus_dVa(:,~ref)));
    dg_dx = dg_dx + sparse(rows,cols,values,nx,nx);
    % dQ_dtheta
    [rows,cols,values] = find(imag(dSbus_dVa(pq,~ref)));
    dg_dx = dg_dx + sparse(rows+nBus,cols,values,nx,nx);
    % dP_dVmag
    [rows,cols,values] = find(real(dSbus_dVm(:,pq)));
    dg_dx = dg_dx + sparse(rows,cols+(nBus-1),values,nx,nx);
    % dQ_dVmag
    [rows,cols,values] = find(imag(dSbus_dVm(pq,pq)));
    dg_dx = dg_dx + sparse(rows+nBus,cols+(nBus-1),values,nx,nx);
    % dP_drho
    dg_dx = dg_dx + sparse(dP_rows(ramp_gen),ix.rho,-PartFact(ramp_gen),nx,nx);
    % fix the derivatives with ZIP[E] contributions
    rows = find(S_E(pq)); cols = find(S_E(pq));
    dg_dx = dg_dx + sparse(rows,cols+(nBus-1),Sd(rows,5).*real(Sd(rows,4)).*Vmag(rows).^(Sd(rows,5)-1),nx,nx);
    dg_dx = dg_dx + sparse(rows+npq,cols+(nBus-1),Sd(rows,5).*imag(Sd(rows,4)).*Vmag(rows).^(Sd(rows,5)-1),nx,nx);
    % fix the derivatives with [Z]IPE contributions
    dg_dx = dg_dx + sparse(rows,cols+(nBus-1),2.*real(Sd(rows,1)).*Vmag(rows),nx,nx);
    dg_dx = dg_dx + sparse(rows+npq,cols+(nBus-1),2.*imag(Sd(rows,1)).*Vmag(rows),nx,nx);
    % fix the derivatives with Z[I]PE contributions
    dg_dx = dg_dx + sparse(rows,cols+(nBus-1),real(Sd(rows,2)),nx,nx);
    dg_dx = dg_dx + sparse(rows+npq,cols+(nBus-1),imag(Sd(rows,2)),nx,nx);
    
    if ~isreal(dg_dx)
        error('dg_dx should be a real matrix.');
    end
end

