clc;clear;


elements = 1:5:100
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
    tic; x = sparsegepp2(K,F,0); t_ge_s(i) = toc;
    err_ge_s(i,1) = ymax-x(id); 
    tic; x = gaussianelimination(K,F); t_ge_f(i) = toc;
    err_ge_f(i,1) = ymax-x(id);
    
    % givens 
    tic; [xf,Af] = Givens_full(K,F,3);  t_givens_f(i) = toc;
    err_givens_f(i,1) = ymax - xf(id);
    tic; [x_COO,A_COO] = Givens_COO(sK,sF); t_givens_s(i) = toc;
    xC = sparse2full(x_COO,'COO');
    err_givens_s(i,1) = ymax - xC(id);
    
    % cholesky
    tic; x = sparseCholeskyFCT(K,F);t_cholesky_s(i) = toc;
    err_cholesky_s(i) = ymax-x(id);
    
    err_cholesky_f(i) = 0;
    
    tic; x = K\F; t_matlab_f(i) = toc;
    err_matlab_f(i) = ymax - x(id); 
    tic; x = sparse(K)\sparse(F); t_matlab_s(i) = toc;
    err_matlab_s(i) = ymax - x(id);

    % get condition number for this stiffness matrix
    condition(i) = cond(K,'fro'); 
    
   
end


figure(1);clf;hold on;
title('error (sparse) vs number of elements');
plot(elements,abs(err_ge_s),'o-');
plot(elements,abs(err_givens_s),'o-');
plot(elements,abs(err_cholesky_s),'o-');
plot(elements,abs(err_matlab_s),'o-');
L=legend('gaussian elimination', 'givens', 'cholesky', 'built in matlab');
L.Location = 'northwest';
ylabel('error');
xlabel('number of elements');

figure(2);clf;hold on;
title('error (full) vs number of elements');
plot(elements,abs(err_ge_f),'o-');
plot(elements,abs(err_givens_f),'o-');
plot(elements,abs(err_cholesky_f),'o-');
plot(elements,abs(err_matlab_f),'o-');
L=legend('gaussian elimination', 'givens', 'cholesky', 'built in matlab');
L.Location = 'northwest';
ylabel('error');
xlabel('number of elements');

figure(3);clf;hold on;
title('error (sparse) vs condition number');
plot(condition,abs(err_ge_s),'o-');
plot(condition,abs(err_givens_s),'o-');
plot(condition,abs(err_cholesky_s),'o-');
plot(condition,abs(err_matlab_s),'o-');
L=legend('gaussian elimination', 'givens', 'cholesky', 'built in matlab');
L.Location = 'northwest';
ylabel('error');
xlabel('condition number');

figure(4);clf;
title('condition number of matrix vs number of elements');
plot(condition,elements);
ylabel('condition number');
xlabel('number of elements');

figure(5);clf; hold on;
title('time vs elements');
plot(elements, t_ge_s, '-o');
plot(elements, t_givens_s, '-o');
plot(elements, t_cholesky_s, '-o');
plot(elements, t_matlab_s, '-o');
ylabel('time (s)');
xlabel('number of elements');
L=legend('gaussian elimination', 'givens', 'cholesky', 'built in matlab');
L.Location = 'northwest';
