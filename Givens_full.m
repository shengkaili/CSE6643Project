function [x,A] = Givens_full(A,b,bw)
%GIVENS_QR solve sparse matrix Ax = b, A & b will be in full matrix format;
% Here A is full ranked
%   1. Apply QR decomposition on A: A = QR;
%   2. solve the triadiagonal matrix: Q'Ax = Q'b
%input:
%   A,b for system Ax = b
%   bw: bw <= 0: give no bandwith information; bw>0: assume certain bandwith,
%   all entries beyond the bandwith will be ingnored.
%output:
%   x: solution of Ax = b;
    [m,n] = size(A);
    % QRD
 if bw <= 0
    for j = 1:n
        i = find(A(:,j),1,'last');
        while  i > j
            if (A(i,j)~=0)
                k = i-1;
                while (k >= j && A(k,j) == 0)
                    k = k-1;
                end
                [c,s] = givens(A(k,j),A(i,j));% 6 flops
                A([k,i],(A(k,:)~=0 |A(i,:)~=0)) = [c,s;-s,c]'*A([k,i],(A(k,:)~=0 |A(i,:)~=0));%
                b([k,i]) = [c,s;-s,c]'*b([k,i]);
                i = k;
            else
                i = i-1;
            end
        end
    end
 else
     for j = 1:n
        mm = min(m,j+bw);
        i = mm;
        while  i > j
            if (A(i,j)~=0)
                k = i-1;
                while (k >= j && A(k,j) == 0)
                    k = k-1;
                end
                [c,s] = givens(A(k,j),A(i,j));% 6 flops
                A([k,i],j:min(n,mm+bw)) = [c,s;-s,c]'*A([k,i],j:min(n,mm+bw));% 6(min(n,mm+bw)-j+1)
                b([k,i]) = [c,s;-s,c]'*b([k,i]);% 6 flops
                i = k;
            else
                i = i-1;
            end
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


