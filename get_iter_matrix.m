function [G] = get_iter_matrix(a, c, b)
% GET_ITER_MATRIX function returns the matrix of iteration (G) in Jacobi
% method for given diagonals (a, b, c) of the tridiagonal matrix of
% a linear equation. Firstly, the function transforms a complex system 
% into a real system according to the formula: 
% C = A+i*B; => M = [A, -B; B, A] (size: 2*Nx2*N).
% Hovever G isn't a typical matrix of iteration. 
% Because matrix of a linear equation is tridiagonal, G has always 
% 8 columns and 2*N rows and stores only these diagonals 
% of a original matrix of iteration which are not zero-diagonals (each 
% diagonal in one column of G). 

% real and imaginary parts 
a_r = real(a); a_i = imag(a);
c_r = real(c); c_i = imag(c);
b_r = real(b); b_i = imag(b);
N = length(c);

% columns of G 
col1 = [0; a_r./c_r(2:N)];
col2 = [b_r./c_r(1:N-1); 0];
col3 = [0; - a_i./c_r(2:N)];
col4 = -c_i./c_r;
col5 = [-b_i./c_r(1:N-1); 0];

G = -[col1, col2, col3, col4, col5];

end

