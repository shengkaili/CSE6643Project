function [x,A] = Givens_full(A,b)
%GIVENS_QR solve sparse matrix Ax = b, A & b will be in full matrix format;
% Here A is full ranked
%   1. Apply QR decomposition on A: A = QR;
%   2. solve the triadiagonal matrix: Q'Ax = Q'b
%input:
%   A,b for system Ax = b
%output:
%   x: solution of Ax = b;
    [m,n] = size(A);
    % QRD
    for j = 1:n
        i = m;
        while  i > j
            if (A(i,j)~=0)
                k = i-1;
                while (k >= j && A(k,j) == 0)
                    k = k-1;
                end
                [c,s] = givens(A(k,j),A(i,j));
                A([k,i],j:n) = [c,s;-s,c]'*A([k,i],j:n);
                b([k,i]) = [c,s;-s,c]'*b([k,i]);
                i = k;
            else
                i = i-1;
            end
        end
    end
    if A(n,n) == 0
        disp('ERROR: rank(A) should equal to n')
        x = NaN;
    else
        % solve x;
        b(n) = b(n)/(A(n,n));
        for i = n-1:-1:1
            b(1:i) = b(1:i) - A(1:i,i+1)*b(i+1);
            b(i) = b(i)/A(i,i);
        end
        x = b(1:n);
    end
end


