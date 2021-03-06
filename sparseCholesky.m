% sparse for both factorization and substitutions
clear all;
nodes = 100; N=2*nodes;
[K,F,ymax,id] = mkbeamproblem(nodes);

% Problem: solve x in Kx=F
% Cholesky decomposition
for j=N:-1:1
    s=0;
    for k=j+1:min(N,j+3)
        s=s+r(j,k)*r(j,k);
    end % 6 flops
    r(j,j)=sqrt(K(j,j)-s);
    
    for i=max(1,j-3):j-1
        s=0;
        for k=j+1:min(N,j+3)
            s=s+r(i,k)*r(j,k);
        end
        r(i,j)=(K(i,j)-s)/r(j,j); 
    end % 8*3=24 flops
end % O(30n flops)

% Now we have K=r*r', so r*r'*x=F
% Use backward substituion to solve r*y=F  (where y=r'*x)
% only nonzero F: F(N-1)=-1;
y(N)=0;
y(N-1)=-1/r(N-1,N-1);
for j=N-2:-1:1
    y(j)=-r(j,j+1:min(N,j+3))*y(j+1:min(N,j+3))'/r(j,j);
end
% O(7n) flops

% Use forward substitution to solve r'*x=y
x(1)=y(1)/r(1,1);
for j=2:N
    x(j)=(y(j)-r(max(1,j-3):j-1,j)'*x(max(1,j-3):j-1)')/r(j,j);
end
% O(7n) flops
% total: O(44n) flops

% visualize the beam
figure,plot((1:2:length(x))/length(x),x(1:2:end),'linewidth',2);
axis equal;
xlabel('x'); ylabel('y'); title('Beam deformation');

errS = ymax - x(id);
xMatlab=K\F;
err = ymax - xMatlab(id);
norm(xMatlab-x') % check the difference b/t sparse Cholesky and in-built solver