function x = gaussianelimination( A,b )

% assume A is square
n = size(A,1);
% initialize gaussian elimination M matrix
I = eye(n,n);
tau = zeros(n,n);

for i = 1:n
   % get value we use to enhilate
   for j = i+1:n
       tau(j,i) = A(j,i)/A(i,i);
   end
   
   % find the multiplier to enhilate the below values
   M = I - tau(:,i)*I(:,i)';
   A = M*A;
   
end

% get the lower matrix
L = I + tau;

y = frdsub(L,b);
x = bcksub(A,y);

function x = bcksub(U,b)
% give me a upper triangular matrix and solution vector b and I can solve
% for input x
%
% assume that U is square

n = size(b,1);
x = zeros(n,1);

x(n) = b(n)/U(n,n);
for i = n-1:-1:1
    x(i) = b(i);
    for j = i+1:n
        x(i) = x(i) - U(i,j)*x(j);
    end
    x(i) = x(i)/U(i,i);
end
function x = frdsub(L,b)
% give me a lower triangular matrix and solution vector b and I can solve
% for input x
%
% assume that L is square

n = size(b,1);
x = zeros(n,1);

x(1) = b(1)/L(1,1);
for i = 2:n
    x(i) = b(i);
    for j = 1:i-1
        x(i) = x(i) - L(i,j)*x(j);
    end
    x(i) = x(i)/L(i,i);
end



