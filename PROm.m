function PROm(~)
%
%

clc;
for t = 1: 30
    maxNumCompThreads(1);
end

% load data, initialize basic variable
Y = importdata('Y.mat');
Phi = {@(z) log(z+1); @(z) atan(z)};
stdGau = [6 6];  % variance of Gaussian
side = 'right';  % left, right, both

% optimizing U, Gamma, V, Lambda
Lst = length(Phi);
N = size(Y(:, :, 1), 2);
Tdpca_Mure = cell(N, 1);
for D2 = 1: N
    [U, Gam, ab, V, Lam, cd] = initPars(Y, [0 1 0], D2);
    M_2dpca = cell(Lst+1, 10);
    M_2dpca{1, 1} = 'Data&STD&SIDE';
    M_2dpca{2, 1} = Y;
    M_2dpca{3, 1} = sprintf('std = [%d, %d], side = %s', stdGau(1), stdGau(2), side);
    M_2dpca{1, 2} = 'phi';
    M_2dpca{1, 3} = 'U';
    M_2dpca{1, 4} = 'Gam';
    M_2dpca{1, 5} = 'ab';
    M_2dpca{1, 6} = 'V';
    M_2dpca{1, 7} = 'Lam';
    M_2dpca{1, 8} = 'cd';
    M_2dpca{1, 9} = 'J';
    M_2dpca{1, 10} = 'Time(s)';
    for lst = 1: length(Phi)
        phi = Phi{lst};
        tic
        [U, Gam, ab, V, Lam, cd, J] = UGVL(Y, stdGau(lst), phi, side, U, Gam, ab, V, Lam, cd);  % left, right, both
        tim = toc;
        
        M_2dpca{lst+1, 2} = sprintf('%s', char(phi));
        M_2dpca{lst+1, 3} = U;
        M_2dpca{lst+1, 4} = Gam;
        M_2dpca{lst+1, 5} = ab;
        M_2dpca{lst+1, 6} = V;
        M_2dpca{lst+1, 7} = Lam;
        M_2dpca{lst+1, 8} = cd;
        M_2dpca{lst+1, 9} = J;
        M_2dpca{lst+1, 10} = tim;
    end
    Tdpca_Mure{D2} = M_2dpca;
    clear M_2dpca
end
save('Tdpca_Mure.mat', 'Tdpca_Mure')
end


