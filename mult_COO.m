function C = mult_COO(A,B)
%MULT_COO calculate the multiplication of matrix in COO format
% C = AB, all matrix are in COO format
NN1 = unique(A(:,1));
NN2 = unique(B(:,2));
if (~isempty(NN1) || ~isempty(NN2))
    C = nan(length(NN1)*length(NN2),3);
    k = 0;
    for j = NN2'
        for i = NN1'
            Ia = find(A(:,1) == i);
            Ib = find(B(:,2) == j);
            Ja = A(Ia,2);
            Jb = B(Ib,1);
            [X,ia,ib] = intersect(Ja,Jb);
            if ~isempty(X)
                k = k+1;
                C(k,1) = i;
                C(k,2) = j;
                aa = A(Ia,3);
                bb = B(Ib,3);
                C(k,3) = aa(ia)'*bb(ib);
            end
        end
    end
    if k == 0
        C = [];
    else
        C = C(1:k,:);
    end
else
    C = [];
    disp('Warning: mult_COO: A or B is empty! Returning [] for C;')
end
end

