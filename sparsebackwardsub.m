function x = sparsebackwardsub(A,b)

n = max(A(:,2));
x = full2sparse(ones(n,1),'COO');

for j = n:-1:1
   Aval = A((A(:,1)==j & A(:,2)==j),3);
   bval = b((
    for i = j:n
     
     
    
   end
end



keyboard

