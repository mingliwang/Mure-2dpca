function Rtnpar = lq_GL(Y, varGau, U, Gam, ab, V, Lam, cd, isLRB, GamOrLam)
%LQ_GL Least Quadratic function for Gam and Lam
%


PREC = 1e-32;
P = size(Y, 3);
switch GamOrLam
    case 'Gam'
        W = eye(size(Y, 2));
        Gam = 0.*Gam;
        for i = 1: P
            Yi = Y(:, :, i);
            
            if isLRB(1)
                Wi = W;
            else
                Wi = V*diag(Lam(:, i))*V';
            end
            UtYi = U'*Yi;
            Ai = diag(0.5*(UtYi*(Wi^2)*UtYi'));
            Bi = -(diag(UtYi*Wi*UtYi') - varGau*trace(Wi));
            aosi = -Bi./(2*Ai+PREC);
            
            gam =  aosi;
            gam(gam<ab(1)) = 0;
            gam(gam>=ab(2)) = 1;
            Gam(:, i) = gam;
        end
        Rtnpar = Gam;
    case 'Lam'
        L = eye(size(Y, 1));
        Lam = 0.*Lam;
        for i = 1: P
            Yi = Y(:, :, i);
            
            if isLRB(2)
                Li = L;
            else
                Li = U*diag(Gam(:, i))*U';
            end
            VtYit = V'*Yi';
            Ai = diag(0.5*(VtYit*(Li^2)*VtYit'));
            Bi = -(diag(VtYit*Li*VtYit') - varGau*trace(Li));
            aosi = -Bi./(2*Ai+PREC);
            
            lam =  aosi;
            lam(lam<cd(1)) = 0;
            lam(lam>=cd(2)) = 1;
            Lam(:, i) = lam;
        end
        Rtnpar = Lam;
end
end

