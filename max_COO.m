function [M, K] = max_COO(A,dim)
%MAX_COO calculate the maximum element of matrix A in COO format
%input:
%   A: matrix for calculation, in COO format
%   dim:dimension to operate along
%output:
%   M: the maximum value(s), in COO format
%   K: row index(es) of maximum value(s) in A. M = A(K,3),
%   the index in the full matrix will be [i,j] = [A(K,1),A(K2)]. K itself
%   will also be in COO format
switch dim
    case 1
        NN = unique(A(:,2));
        M = nan(length(NN),3);
        K = nan(length(NN),3);
        M(:,1) = 1;
        K(:,1) = 1;
        k = 0;
        for i = NN'
            k = k+1;
            M(k,2) = i;
            K(k,2) = i;
            II = find(A(:,2) == i);
            aa = A(II,3);
            [mm, ii] = max(aa);
            M(k,3) = mm;
            K(k,3) = II(ii);
        end
    case 2
        NN = unique(A(:,1));
        M = nan(length(NN),3);
        K = nan(length(NN),3);
        M(:,2) = 1;
        K(:,2) = 1;
        k = 0;
        for i = NN'
            k = k+1;
            M(k,1) = i;
            K(k,1) = i;
            II = find(A(:,1) == i);
            aa = A(II,3);
            [mm, ii] = max(aa);
            M(k,3) = mm;
            K(k,3) =II(ii);
        end
    otherwise
        disp('ERROR: max_COO: dim cannot be other than 1 or 2')
        keyboard
end
end

