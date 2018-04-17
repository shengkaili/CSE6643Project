%% test the givens solvers for the problem
clear all;
nodes = 500; N=2*nodes;
[K,F,ymax,id] = mkbeamproblem(nodes);
disp('Applying matlab inverse method ...')
tic;xMatlab=K\F;tM = toc
err = ymax - xMatlab(id);
disp('Applying Givens_full method ...')
tic;[xf,Af] = Givens_full(K,F,3);tf = toc
norm(xMatlab-xf) % check the difference Givens_full and in-built solver
errG_f = ymax - xf(id);
errG_f2 = ymax - xf2(id);
K_COO = full2sparse(K,'COO');
F_COO = full2sparse(F,'COO');
disp('Applying Givens_COO method ...')
tic;[x_COO,A_COO] = Givens_COO(K_COO,F_COO);tC = toc
xC = sparse2full(x_COO,'COO');
AC = sparse2full(A_COO,'COO');
errG_C = ymax - xC(id);
norm(xMatlab-xC) % check the difference Givens_COO and in-built solver