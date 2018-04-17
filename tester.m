clc;clear;


elements = 1:20;
error = zeros(length(elements),1);
condition = error;
for i = 1:length(elements)
   
    % get number of nodes
    n = elements(i);
    % get sample problem
    [K,F,ymax,id] = mkbeamproblem(n);
    % convert K, F into COO sparse format
    sK = full2sparse(K,'COO');
    sF = full2sparse(F,'COO');    
   
    % LU decomposition via guassian elimination
    x = sparsegepp2(K,F,0);
    err_ge_s(i,1) = ymax-x(id);
    x = gaussianelimination(K,F);
    err_ge_f(i,1) = ymax-x(id);
    
    % givens 
    [xf,Af] = Givens_full(K,F,3);
    err_givens_f(i,1) = ymax - xf(id);
    [x_COO,A_COO] = Givens_COO(sK,sF);
    xC = sparse2full(x_COO,'COO');
    err_givens_s(i,1) = ymax - xC(id);
    
    % cholesky
    x = sparseCholeskyFCT(K,F);
    err_cholesky_s = ymax-x(id);

    % get condition number for this stiffness matrix
    condition(i) = cond(K,'fro'); 
    
    
    

    
    
    
    

    
    
    
end

keyboard
