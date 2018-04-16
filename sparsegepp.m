function [A,b] = sparsegepp(A,b)
%SPARSEGEPP the user inputs a square A matrix and solution vector b in
%'COO' sparse format. The output is an upper triangular matrix and properly
%rotated solution vector in 'COO' sparse format. This algorithm is doing
%row pivoting and gaussian elimination to reduce A to upper triangular.

% number of columns in square matrix A
n = max(A(:,2));

% initalize anihilation matrix
M = A;
% set values in M to be NaN to make sure we are eliminating everything we
% expect to
M(:,3) = NaN;

% loop through each column
for i = 1:n
   
    % loop back to create a new subset if we pivoted
     redo = true;
    while redo
        
        % get subset to eliminate, only get items at or below the diagonal
        subsetloc = find(A(:,2)==i & A(:,1)>=A(:,2)); 
        subset = A(subsetloc,:);
        
        % get indexes to compare to pivot
        [~,maxindx] = max(abs(subset(:,3)),[],1);
        maxindx = subset(maxindx,1);
        % pivot if there is a value larger in magnitude than diagonal
        if maxindx ~= i
            temp = A(:,1)==maxindx;
            temp2 = b(:,1)==maxindx;
            b(b(:,1)==i,1) = maxindx;
            M(A(:,1)==i,1) = maxindx;
            A(A(:,1)==i,1) = maxindx;
            b(temp2,1) = i;
            M(temp,1) = i;
            A(temp,1) = i;
            redo = true;
        else
            redo = false;
        end
    end
    
    % diagonal term is always first
    [~,order] = sort(subset(:,1));
    subset = subset(order,:);
    subsetloc = subsetloc(order);
    
    % find our anhilation values and update the A, b matrix
    for j = 2:size(subset,1)
        M(subsetloc(j),3) = -subset(j,3)/subset(1,3);
        
        row = subset(j,1);
        col = subset(j,2);
        
        % update A with M
        for k = i:n
            
            % find index we are updating
            indx = find(A(:,1)==row & A(:,2)==k);
            indx2 = find(A(:,1)==col & A(:,2)==k);
            
            if isempty(indx) && ~isempty(indx2)
                % create sparse A entry
                A(end+1,:) = [row, k, 0];
                indx = size(A,1);
            end
            if ~isempty(indx2)
                A(indx,3) = A(indx,3) + M(subsetloc(j),3)*A(indx2,3);
            end
        end
        
        % update F with M
        indx3 = find(b(:,1)==row);
        indx4 = find(b(:,1)==col);
        
        if isempty(indx3) && ~isempty(indx4)
            b(end+1,:) = [row, 1, 0];
            indx3 = size(b,1);
        end
        if ~isempty(indx4)
            b(indx3,3) = b(indx3,3) + M(subsetloc(j),3)*b(indx4,3);
        end

        
    end
    
    
    try
    M(subsetloc(1),3) = 1;
    catch
       keyboard 
    end
   

end



