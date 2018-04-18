clc;clear;


elements = 2:10:300;
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
    xnorm_ge_s(i) = norm(x,2);
    x_ge_s(i,1) =x(id);
    err_ge_s(i,1) = ymax-x(id); 
    tic; x = gaussianelimination(K,F); t_ge_f(i) = toc;
    x_ge_f(i,1) =x(id);
    err_ge_f(i,1) = ymax-x(id);
    
    % givens 
    tic; [x,Af] = Givens_full(K,F,3);  t_givens_f(i) = toc;
    x_givens_f(i,1)= x(id);
    err_givens_f(i,1) = ymax - x(id);
    tic; [x,A_COO] = Givens_COO(sK,sF); t_givens_s(i) = toc;
    x = sparse2full(x,'COO');
    x_givens_s(i,1) =x(id);
    xnorm_givens_s(i,1) = norm(x,2);
    err_givens_s(i,1) = ymax - x(id);
    
    % cholesky
    tic; x = sparseCholeskyFCT(K,F);t_cholesky_s(i) = toc;
    x_cholesky_s(i,1)= x(id);
    xnorm_cholesky_s(i,1) = norm(x,2);
    err_cholesky_s(i) = ymax-x(id);
    
    err_cholesky_f(i) = 0;
    
    tic; x = K\F; t_matlab_f(i) = toc;
    x_matlab_f(i,1) = x(id);
    err_matlab_f(i) = ymax - x(id); 
    tic; x = sparse(K)\sparse(F); t_matlab_s(i) = toc; x= full(x);
    xnorm_matlab_s = norm(x,2);
    x_matlab_s(i,1)= x(id);
    err_matlab_s(i) = ymax - x(id);
    
    
    % get condition number for this stiffness matrix
    condition(i) = cond(K,'fro'); 
    
    [~,S,~] = svd(K);
    sigma1(i) = S(end);
    fnorm(i) = norm(F,2);
    
    numnnz(i) = nnz(K);
    numful(i) = size(K,1)*size(K,2)-nnz(K);
    
   
end
condition2 = elements.^(2/3);


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
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')


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
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

figure(3);clf;hold on;
title('error (sparse) vs condition number');
plot(condition2,abs(err_ge_s),'o-');
plot(condition2,abs(err_givens_s),'o-');
plot(condition2,abs(err_cholesky_s),'o-');
plot(condition2,abs(err_matlab_s),'o-');
L=legend('gaussian elimination', 'givens', 'cholesky', 'built in matlab');
L.Location = 'northwest';
ylabel('error');
xlabel('condition number');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

figure(4);clf;hold on;
plot(elements,condition);
plot(elements,condition2);
title('condition number of matrix vs number of elements');
ylabel('condition number');
xlabel('number of elements');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

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
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

% % h
% figure(6);clf; hold on;
% plot(err_ge_s,fnorm./xnorm_ge_s./sigma1);
% plot(err_matlab_s,fnorm./xnorm_matlab_s./sigma1);
% plot(err_givens_s',fnorm./xnorm_givens_s'./sigma1);
% plot(err_cholesky_s,fnorm./xnorm_cholesky_s'./sigma1);
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% legend('condition number', 'plesha condition number');

