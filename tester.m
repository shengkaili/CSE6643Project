
nodes = [15,100];
error = zeros(length(nodes),1);
condition = error;


for i = 1:length(nodes)
   
    % get number of nodes
    n = nodes(i);
    % get sample problem
    [K,F,ymax,id] = mkbeamproblem(n);
    % get condition number for this stiffness matrix
    condition(i) = cond(K,'fro'); 
    
    % convert K, F into COO sparse format
    sK = full2sparse(K,'COO');
    sF = full2sparse(F,'COO');
    
    % LU decomposition via guassian elimination
    [A,b] = sparsegepp(sK,sF);
    
%     % solve sparse matrix by cheating, write backward subsitution
%     x = sparse2full(A,'COO')\sparse2full(b,'COO');
%     %     x = sparsebackwardsub(A,b);
%     
%     error(i) = abs(ymax - x(id))/ymax;
    
    
    

    
    
    
end
