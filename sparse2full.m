function B = sparse2full(A,style)
% as a check, convert the sparse matrix back to full and make sure our
% sparse method is working
%
% this works well for the K matrix but for the force vector we have some
% sparse entries at the boundaries so it doesn't return the full size
% expected

if strcmp(style, 'COO')
    B = zeros(max(A(:,1)), max(A(:,2)));
    for i = 1:size(A,1)
       B(A(i,1),A(i,2)) = A(i,3); 
    end
else
   keyboard 
end