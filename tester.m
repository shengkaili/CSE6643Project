
nodes = [2,5,15,100];


for n = nodes
   
    [K,F,ymax,id] = mkbeamproblem(n);
    
    kappa = cond(K,'fro');
    
    sK = full2sparse(K,'COO');
    sF = full2sparse(F,'COO');
    
    [A,b] = sparsegepp(sK,sF);
    
    x = sparse2full(A,'COO')\sparse2full(b,'COO');
    
%     x = sparsebackwardsub(A,b);
    
    
    
end
