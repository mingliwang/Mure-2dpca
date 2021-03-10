function [U, Gam, ab, V, Lam, cd] = initPars(Y, isLRB, D2)
%INITPARS Initial U, Gam, V, Lam; a, b c, d
%

[M, N, P] = size(Y);
ratio = 1;

if (isLRB(1) || isLRB(3))
    Y2 = reshape(Y, M, [], 1);
    [U, ~, ~] = svd(Y2*Y2');
    
    D1 = min(round(M*ratio), M);  % D1<=M
    U = U(:, 1:D1);
    Gam = ones(D1, P);
    ab = [0.5 0.5];  % 0<=a<=0.5<=b<=1
else
    U = [];
    Gam = [];
    ab = [];
end

if (isLRB(2) || isLRB(3))
    Y3 = reshape(permute(Y, [2 1 3]), N, [], 1)';
    [V, ~, ~] = svd(Y3'*Y3);
    % D2 = min(round(N*ratio), N);  % D2<=N
    V = V(:, 1:D2);
    Lam = ones(D2, P);
    cd = [0.5 0.5];  % 0<=a<=0.5<=b<=1
else
    V = [];
    Lam = [];
    cd = [];
end
end
