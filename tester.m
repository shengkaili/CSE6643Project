
nodes = [3,100,200];
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
    
    
    x = sparsegepp2(K,F,0);
    errorspr(i) = abs(ymax - x(id))/ymax;
    
    x = gaussianelimination(K,F);
    errornrm(i) = abs(ymax - x(id))/ymax;
    
%     % solve sparse matrix by cheating, write backward subsitution
%     x = sparse2full(A,'COO')\sparse2full(b,'COO');
    %     x = sparsebackwardsub(A,b);
%     
    
    
    
    

    
    
    
end

keyboard
