function [x,A] = Givens_COO(A,b)
%GIVENS_QR solve sparse matrix Ax = b, A & b will be in COO format;
%   1. Apply QR decomposition on A: A = QR;
%   2. solve the triadiagonal matrix: Q'Ax = Q'b
%input:
%   A,b for system Ax = b, A & b are saved as the 'COO' style
%output:
%   x: solution of Ax = b, as 'COO' style;
m = max(A(:,1));
n = max(A(:,2));
% QRD
for j = 1:n
    Ka = find(A(:,2) == j);
    Ia = A(Ka,1);
    aa = A(Ka,3);
    [IIa,KKa] = sort(Ia,'descend');
    if isempty(IIa)
        disp('ERROR: Givens_COO: Matrix is not full ranked!')
        keyboard
    end
    i = 1;
    while  i<length(IIa) && IIa(i)>j
        I2 = IIa(i);
        I1 = IIa(i+1);
        [c,s] = givens(aa(KKa(i+1)),aa(KKa(i)));
        aa(KKa(i+1)) = c*aa(KKa(i+1))-s*aa(KKa(i));
        A(Ka(KKa(i+1)),3) = aa(KKa(i+1));
        A(Ka(KKa(i)),3) = 0;
        K1 = find(A(:,1) == I1 & A(:,2) > j);
        K2 = find(A(:,1) == I2 &  A(:,2) > j);
        JJ = unique([A(K1,2);A(K2,2)]);
        for jj = JJ'
            KK1 = find(A(K1,2) == jj);
            if isempty(KK1)
                x1 = 0;
                A = [A;[I1,jj,NaN]];
                kk1 = size(A,1);
            else
                x1 = A(K1(KK1),3);
                kk1 = K1(KK1);
            end
            KK2 = find(A(K2,2) == jj);
            if isempty(KK2)
                x2 = 0;
                A = [A;[I2,jj,NaN]];
                kk2 = size(A,1);
            else
                x2 = A(K2(KK2),3);
                kk2 = K2(KK2);
            end
            A(kk1,3) = c*x1-s*x2;
            A(kk2,3) = s*x1+c*x2;
        end
        if ~isempty(find(b(:,1)==I1,1)) || ~isempty(find(b(:,1)==I2,1))
            kb1 = find(b(:,1)==I1);
            kb2 = find(b(:,1)==I2);
            if isempty(kb1)
                b1 = 0;
                b = [b;[I1,1,NaN]];
                kb1 = size(b,1);
            else
                b1 = b(kb1,3);
            end
            if isempty(kb2)
                b2 = 0;
                b = [b;[I2,1,NaN]];
                kb2 = size(b,1);
            else
                b2 = b(kb2,3);
            end
            b(kb1,3) = c*b1-s*b2;
            b(kb2,3) = s*b1+c*b2;
        end
        i = i+1;
    end
    if i == length(IIa) && IIa(i)<j
        kk = find(A(:,1)==IIa(i));
        A(A(:,1)==j,1) = IIa(i);
        A(kk,1) = j;
    end
end
A(A(:,3)==0,:) = [];
b(A(:,3)==0,:) = [];
% A_f = sparse2full(A,'COO');
% solve x
Ib = find(b(:,1)==n);
Ia = find(A(:,1)==n & A(:,2)==n);
if ~isempty(Ib) && ~isempty(Ia)
    b(Ib,3) = b(Ib,3)/(A(Ia,3));
    for i = n-1:-1:1
        b1 = b((b(:,1)==(i+1)),3);
        if ~isempty(b1)
            for j = 1:i
                A1 = A(A(:,1)==j & A(:,2)==(i+1),3);
                if ~isempty(A1)
                    Ib = find(b(:,1)==j);
                    if isempty(Ib)
                        b = [b;[j,1,0]];
                        Ib = size(b,1);
                    end
                    b(Ib,3) = b(Ib,3) - A1*b1;
                end
            end
        end
        Ib = find(b(:,1) == i);
        if ~isempty(Ib)
           b(Ib,3) = b(Ib,3)/A(A(:,1)==i & A(:,2)==i,3);
        end
    end
end
x = b(b(:,1)<=n,:);
x(x(:,3)==0,:) = [];
end

