Af = rand(5,3);
Bf = rand(3,6);
% Af(2,:) = 0;
% Bf(:,3) = 0;
A = full2sparse(Af,'COO');
B = full2sparse(Bf,'COO');
Cf = Af*Bf;
C = mult_COO(A,B);
err = C-full2sparse(Cf,'COO');
errf = sparse2full(C,'COO') - Cf;

dim = 2;
Af(Af == 0) = NaN;
[Mf,If] = max(Af,[],dim);
[M,K] = max_COO(A,dim);
I = A(K(:,3),dim);

[mf,IIf] = min(Af,[],dim);
[m,KK] = min_COO(A,dim);
II = A(KK(:,3),dim);