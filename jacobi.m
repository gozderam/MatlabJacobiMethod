function [x, iter] = jacobi(a, c, b, d, tol)
% JACOBI function solves M*x = d with Jacobi method, where:
% M - tridiagonal matrix of complex numbers (NxN)
% a, c, b - non-zero diagonals of M
% d - vertical vector (Nx1)
% tol - tolerance 
% function returns: 
% x - solution vector
% iter - number of iterations (if Jacobi method is divergent, function
% shows information)

% d vector and transforming complex vector into real vector
c_r = real(c); 
d = [real(d); imag(d)];
d = d./[c_r; c_r];
N = length(c);

% matrix of iteration - version for tridiagonal matrix
% G1 is used in 1:N elements of x, in order to use it in N+1:2*N elements
% of x we have to move columns and add minus sign
G1 = get_iter_matrix(a, c, b);
G2 = [-G1(:, 3:5), G1(:, 1:2)];
% using to parts of G (to avoid unnecessary multiplications)

% number of iterations
it = 0; 

% initial values
xk = zeros(2*N, 1);
xk1 = ones(2*N, 1);

% jacobi algorithm
while(norm(xk-xk1)>=tol)
     xk = xk1;
    it = it+1;
    % firs row
        xk1(1) = G1(1, :) * [0, xk(2), 0, xk(N+1), xk(N+2)]' + d(1);
    % rows from 2 to N-1
        for i=2:N-1
            x_m=[xk(i-1), xk(i+1), xk(N+i-1), xk(N+i), xk(N+i+1)]';
            xk1(i) = G1(i, :) * x_m + d(i);
        end
    % row N
        xk1(N) = G1(N, :) * [xk(N-1), 0, xk(2*N-1), xk(2*N), 0]' + d(N); 
    % row N+1 
        xk1(N+1) = G2(1, :) * [0, xk(1), xk(2), 0, xk(N+2)]' + d(N+1);
    % rows from N+2 to 2*N-1
        for i=(N+2):(2*N-1)
            x_m = [xk(i-N-1), xk(i-N), xk(i-N+1), xk(i-1), xk(i+1)]';
            xk1(i) = G2(i-N, :) * x_m + d(i);
        end
    % row 2*N
        xk1(2*N) = G2(N, :) * [xk(N-1), xk(N), 0, xk(2*N -1), 0]' + d(2*N);
      
    % if the method is divergent
     if(it>=100000 || norm(xk1) > eps^(-1/2))
        iter = it;
        x = xk1(1:N) + 1i*xk1(N+1:2*N);
        disp("Function stopped - in this case Jacobi method is divergent.");
        return;
     end
        
end

% returned values
x = xk1(1:N) + 1i*xk1(N+1:2*N);
iter = it;
end


