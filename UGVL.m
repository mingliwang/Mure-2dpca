function [U, Gam, ab, V, Lam, cd, J] = UGVL(Y, stdGau, phi, side, U, Gam, ab, V, Lam, cd)
%UGVL Updating U, Gamma, V, Lambda in succession
%

isLRB = [~isempty(regexpi(side, 'left')), ~isempty(regexpi(side, 'right')), ...
    ~isempty(regexpi(side, 'both'))];

varGau = stdGau^2;
phid = fhd(phi);
J0 = evalOF(Y, varGau, phi, U, Gam, V, Lam, isLRB);

LOOP = 60;
J = zeros(LOOP, 1);
for loop = 1: LOOP
    if (isLRB(1) || isLRB(3))
        U = pd_UV(Y, varGau, phid, U, Gam, V, Lam, isLRB, 'U');
        Gam = lq_GL(Y, varGau, U, Gam, ab, V, Lam, cd, isLRB, 'Gam');
    end
    if (isLRB(2) || isLRB(3))
        V = pd_UV(Y, varGau, phid, U, Gam, V, Lam, isLRB, 'V');
        Lam = lq_GL(Y, varGau, U, Gam, ab, V, Lam, cd, isLRB, 'Lam');
    end
    J(loop) = evalOF(Y, varGau, phi, U, Gam, V, Lam, isLRB);
    if (loop>=2) && (abs(J(loop)-J(loop-1))<=1e-9)
        J = J(1:loop);
        break
    end
end
J = cat(1, J0, J);
end


function J = evalOF(Y, varGau, phi, U, Gam, V, Lam, isLRB)
%EVALOF Evaluate objective function
%

L = eye(size(Y, 1));
W = eye(size(Y, 2));

P = size(Y, 3);
phiE = zeros(P, 1);
for i = 1: P
    Yi = Y(:, :, i);
    
    if isLRB(1)
        Li = U*diag(Gam(:, i))*U';
        Wi = W;
    elseif isLRB(2)
        Wi = V*diag(Lam(:, i))*V';
        Li = L;
    elseif isLRB(3)
        Li = U*diag(Gam(:, i))*U';
        Wi = V*diag(Lam(:, i))*V';
    end
    Ei = 0.5*sum(sum((Yi - Li*Yi*Wi).^2)) + varGau*trace(Li)*trace(Wi);
    phiE(i) = phi(Ei);
end
J = sum(phiE);
end


function phid = fhd(phi)
%FHD derivative of function handle
%

syms z
phid = eval(['@(z)' vectorize(char(diff(phi(z))))]);
% or using: phid = matlabFunction(diff(phi(z)));
end


