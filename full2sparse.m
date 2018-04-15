function A = full2sparse(A,style)
% input a full matrix or vector "A" and outputs a sparse version in "COO"
% format (row, column, value)
% 
% https://en.wikipedia.org/wiki/Sparse_matrix
%
% this format isn't set in store, just wanted to implement something
%

if strcmp(style, 'COO')
    indx = find(A~=0);
    [row, col] = ind2sub(size(A),indx);
    A = [row, col, A(indx)];
else
    keyboard
end


