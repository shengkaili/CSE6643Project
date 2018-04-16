%% test the givens solvers for the problem
clear all;
nodes = 100; N=2*nodes;
[K,F,ymax,id] = mkbeamproblem(nodes);
disp('Applying Givens_full method ...')
tic;[xf,Af] = Givens_full(K,F);toc
errG_f = ymax - xf(id);
K_COO = full2sparse(K,'COO');
F_COO = full2sparse(F,'COO');
disp('Applying Givens_COO method ...')
tic;[x_COO,A_COO] = Givens_COO(K_COO,F_COO);toc
xC = sparse2full(x_COO,'COO');
AC = sparse2full(A_COO,'COO');
errG_C = ymax - xC(id);
disp('Applying matlab inverse method ...')
tic;xMatlab=K\F;toc
err = ymax - xMatlab(id);
norm(xMatlab-xf) % check the difference Givens_full and in-built solver
norm(xMatlab-xC) % check the difference Givens_COO and in-built solver