function x = Givens_tridiag(A,b,n)
%GIVENS_TRIDIAG solve sparse symetric matrix:
%   1. transpform the matrix into tridiagonal matrix: A = Q'TQ;
%   2. solve the triadiagonal matrix: Dy = Qb, where y = Qx;
%   3. x = Q'y.
%input:
%   A,b for system Ax = b, A is symetric and saved as the 'COO' style
%   n: n = 0: solve the whole matrix; n>0: assume the bandwith of
%   matrix A is n.
%output:
%   x: solution of Ax = b;
Nmax = max(A(:,1));

for i = 1:Nmax
    


end
end


