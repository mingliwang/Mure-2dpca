function Rtnpar = pd_UV(Y, varGau, phid, U, Gam, V, Lam, isLRB, UorV)
%PD_UV Polor Decomposition for U and V.
%

P = size(Y, 3);
mxIter = 1e3;
cvgPrc = 1e-5;
switch UorV
    case 'U'
        Uupd = PD_U(Y, varGau, phid, U, Gam, V, Lam, isLRB, P, mxIter, cvgPrc);
        Rtnpar = Uupd;
    case 'V'
        Vupd = PD_V(Y, varGau, phid, U, Gam, V, Lam, isLRB, P, mxIter, cvgPrc);
        Rtnpar = Vupd;
end
end


function Uupd = PD_U(Y, varGau, phid, U, Gam, V, Lam, isLRB, P, mxIter, cvgPrc)
%PD_U
%


W = eye(size(Y, 2));
bsGam = logical(sum(Gam, 2))';
for iter = 1: mxIter
    nblUJ = 0;
    for i = 1: P
        Yi = Y(:, :, i);
        
        Li = U*diag(Gam(:, i))*U';
        if isLRB(1)
            Wi = W;
        else
            Wi = V*diag(Lam(:, i))*V';
        end
        Ei = 0.5*sum(sum((Yi - Li*Yi*Wi).^2)) + varGau*trace(Li)*trace(Wi);
        nblEiphi = phid(Ei);
        nblUEi = (Yi*(Wi^2)*Yi')*(U*diag(Gam(:, i).^2)) - 2*(Yi*Wi*Yi')*(U*diag(Gam(:, i)));
        
        nblUJ = nblUJ + nblEiphi*nblUEi;
    end
    nblUJ = nblUJ/(mean(abs(nblUJ(:))) + 1e-6);
    [tpV, tpS] = eig(nblUJ'*nblUJ);  % single precision may produce negative eigenvalues
    tpS_isr = invSR(tpS);
    Uupd = -(nblUJ*(tpV*tpS_isr*tpV'));
    % [tpU, ~, tpV] = svd(nblUJ, 'econ');
    % Uupd = -tpU*tpV';
    
    if cvgCnd(Uupd, U, bsGam, cvgPrc)
        break
    else
        U = Uupd;
    end
end
end


function Vupd = PD_V(Y, varGau, phid, U, Gam, V, Lam, isLRB, P, mxIter, cvgPrc)
%PU_V iterating updata V
%   Analogous to PD_U

L = eye(size(Y, 1));
bsLam = (logical(sum(Lam, 2)))';
for iter = 1: mxIter
    nblVJ = 0;
    for i = 1: P
        Yi = Y(:, :, i);
        
        if isLRB(2)
            Li = L;
        else
            Li = U*diag(Gam(:, i))*U';
        end
        Wi = V*diag(Lam(:, i))*V';
        Ei = 0.5*sum(sum((Yi - Li*Yi*Wi).^2)) + varGau*trace(Li)*trace(Wi);
        nblEiphi = phid(Ei);
        nblVEi = (Yi'*(Li^2)*Yi)*(V*diag(Lam(:, i).^2)) - 2*(Yi'*Li*Yi)*(V*diag(Lam(:, i)));
        
        nblVJ = nblVJ + nblEiphi*nblVEi;
    end
    nblVJ = nblVJ/(mean(abs(nblVJ(:))) + 1e-6);
    [tpV, tpS] = eig(nblVJ'*nblVJ);
    tpS_isr = invSR(tpS);
    Vupd = -(nblVJ*(tpV*tpS_isr*tpV'));
    
    if cvgCnd(Vupd, V, bsLam, cvgPrc)
        break
    else
        V = Vupd;
    end
end
end


function tpS_isr = invSR(tpS)
%INVSR Inverse square root of tpS
%

s = diag(tpS);
ind = (s>0);
s(s<=0) = 1e16;
s = sqrt((s.^(-1)));
s = s.*ind;
tpS_isr = diag(s);
end


function isCvg = cvgCnd(New, Old, colVld, cvgPrc)
%CVGCOND Convergence condition
%

error = sum(abs((sum(New.*Old, 1) - sum(Old.^2, 1)).*colVld));
isCvg = (error <= cvgPrc);
end

