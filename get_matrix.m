function [M, d] = get_matrix(N, spectral_radius_greater_than_1)
% GET_MATRIX - for test only. M - random tridiagonal matrix of a linear 
% equation Mx=d. N - size of M.
% spectral_radius_greater_than_1 == 1 =>  spectral radius of M is greater
% than 1 => Jacobi method is divergent
% spectral_radius_greater_than_1 == 0 =>  spectral radius of M is less
% than 1 => Jacobi method is covergent

d = rand(N, 1) + 1i*rand(N, 1);
a = rand(N-1, 1) + 1i*rand(N-1, 1);
c = rand(N, 1) + 1i*rand(N, 1);
b = rand(N-1, 1) + 1i*rand(N-1, 1);
M = diag(a, -1) + diag(c) + diag(b, 1);
M_jacobi = [real(M), -imag(M); imag(M), real(M)];
D = diag(diag(M_jacobi));
D = inv(D);
M_jacobi = -D * ( M_jacobi - diag(diag(M_jacobi)));

if(spectral_radius_greater_than_1 == 0)
    while (max(abs(eig(M_jacobi))) >=1 )
    a = rand(N-1, 1) + 1i*rand(N-1, 1);
    c = rand(N, 1) + 1i*rand(N, 1);
    b = rand(N-1, 1) + 1i*rand(N-1, 1);
    M = diag(a, -1) + diag(c) + diag(b, 1);
    M_jacobi = [real(M), -imag(M); imag(M), real(M)];
    D = diag(diag(M_jacobi));
    D = inv(D);
    M_jacobi = -D * ( M_jacobi - diag(diag(M_jacobi)));
    end
elseif(spectral_radius_greater_than_1 == 1)
    while (max(abs(eig(M_jacobi))) <=1 )
    a = rand(N-1, 1) + 1i*rand(N-1, 1);
    c = rand(N, 1) + 1i*rand(N, 1);
    b = rand(N-1, 1) + 1i*rand(N-1, 1);
    M = diag(a, -1) + diag(c) + diag(b, 1);
    M_jacobi = [real(M), -imag(M); imag(M), real(M)];
    D = diag(diag(M_jacobi));
    D = inv(D);
    M_jacobi = -D * ( M_jacobi - diag(diag(M_jacobi)));
    end
end  
end

