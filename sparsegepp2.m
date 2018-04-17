function x = sparsegepp2(A,b,okpivot)

% number of columns in square matrix A
n = size(A,1);
L = eye(n,n);
U = eye(n,n);

for j = 1:n
    
    % solve Ljuj=aj ->uj
    % 1. determine positions of nonzero uj and order to solve
    % 2. solve    
    Lj = L(1:j-1,1:j-1);
    aj = A(1:j-1,j);
    uj = solveLx(Lj,aj);
    U(1:j-1,j) = uj;
    
    Ljp = L(j:n,1:j-1);
    ajp = A(j:n,j);
    bj = ajp - Ljp*uj;
    
    %     pivot
    if okpivot
        [m,i] = max(abs(bj));
        if i+j-1~=j
            temp = bj(1);
            bj(1) = m;
            bj(i) = temp;
            
            P = eye(n,n);
            temp = P(j,:);
            P(j,:) = P(i+j-1,:);
            P(i+j-1,:) = temp;
            
            A = P*A;
            b = P*b;
            
        end
    end
    ujj = bj(1);
    U(j,j) = ujj;
    L(j:n,j) = bj/ujj;  
end 

y = L\b;
x = U\y;

function x = solveLx(L,b)
x = b;
n = length(x);
for j = 1:n
    if (x(j)~=0)
        x(j+1:n) = x(j+1:n) - L(j+1:n,j)*x(j);
    end
end

function x = solveUx(U,b)
x = b;
n = size(x,1);
for j = n:-1:1
        x(j) = (x(j) - U(j,j+1:n)*x(j+1:n))/U(j,j);
end









