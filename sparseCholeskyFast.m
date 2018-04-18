function x= sparseCholeskyFast(K,F)
% no max or min functions
% Problem: solve x in Kx=F
% Cholesky decomposition
N = length(F);
% j=N
if N<=3
    x= sparseCholeskyFCT(K,F);
else
    r(N,N)=sqrt(K(N,N));
    for i=(N-3):(N-1)
        r(i,N)=(K(i,N))/r(N,N);
    end % 8*3=24 flops
    
    % j=N-1
    s=r(N-1,N)*r(N-1,N);
    r(N-1,N-1)=sqrt(K(N-1,N-1)-s);
    for i=(N-4):(N-2)
        s=r(i,N)*r(N-1,N);
        r(i,N-1)=(K(i,N-1)-s)/r(N-1,N-1);
    end % 8*3=24 flops
    
    % j=N-2
    s=0;
    for k=N-1:N
        s=s+r(N-2,k)*r(N-2,k);
    end % 6 flops
    r(N-2,N-2)=sqrt(K(N-2,N-2)-s);
    for i=(N-5):(N-3)
        s=0;
        for k=N-1:N
            s=s+r(i,k)*r(N-2,k);
        end
        r(i,N-2)=(K(i,N-2)-s)/r(N-2,N-2);
    end % 8*3=24 flops
    
    
    for j=N-3:-1:4
        s=0;
        for k=(j+1):(j+3)
            s=s+r(j,k)*r(j,k);
        end % 6 flops
        r(j,j)=sqrt(K(j,j)-s);
        
        for i=(j-3):(j-1)
            s=0;
            for k=(j+1):(j+3)
                s=s+r(i,k)*r(j,k);
            end
            r(i,j)=(K(i,j)-s)/r(j,j);
        end % 8*3=24 flops
    end % O(30n flops)
    
    % j=3
    s=0;
    for k=4:6
        s=s+r(3,k)*r(3,k);
    end % 6 flops
    r(3,3)=sqrt(K(3,3)-s);
    for i=1:2
        s=0;
        for k=4:6
            s=s+r(i,k)*r(3,k);
        end
        r(i,3)=(K(i,3)-s)/r(3,3);
    end % 8*3=24 flops
    
    % j=2
    s=0;
    for k=3:5
        s=s+r(2,k)*r(2,k);
    end % 6 flops
    r(2,2)=sqrt(K(2,2)-s);
    s=0;
    for k=3:5
        s=s+r(1,k)*r(2,k);
    end
    r(1,2)=(K(1,2)-s)/r(2,2);
    
    % j=1
    s=0;
    for k=2:4
        s=s+r(1,k)*r(1,k);
    end % 6 flops
    r(1,1)=sqrt(K(1,1)-s);
    
    % Now we have K=r*r', so r*r'*x=F
    % Use backward substituion to solve r*y=F  (where y=r'*x)
    % only nonzero F: F(N-1)=-1;
    y(N)=0;
    y(N-1)=-1/r(N-1,N-1);
    y(N-2)=-r(N-2,N-1:N)*y(N-1:N)'/r(N-2,N-2);
    for j=N-3:-1:1
        y(j)=-r(j,(j+1):(j+3))*y((j+1):(j+3))'/r(j,j);
    end
    % O(7n) flops
    
    % Use forward substitution to solve r'*x=y
    x(1)=y(1)/r(1,1);
    x(2)=(y(2)-r(1,2)*x(1))/r(2,2);
    x(3)=(y(3)-r(1:2,3)'*x(1:2)')/r(3,3);
    for j=4:N
        x(j)=(y(j)-r((j-3):(j-1),j)'*x((j-3):(j-1))')/r(j,j);
    end
end
end

