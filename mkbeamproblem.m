function [K,F,ymax,id] = mkbeamproblem(n)
% MKBEAMPROBLEM(nodes), input the number of nodes you want the FEM to have.
% This will output the stiffness matrix K, force vector F, analytical end
% position of the beam ymax and location in solution vector x (K\F) that
% should be compared to the analytical solution. 

% defining variables..
s = 0.125; % 1 square inch section
L = 1; % length (m)
E = 68.9E9; % elasticity of AL6061 (Pa)
I = (s*0.0254)^4/12; % area moment of inertia (m^4)
P = 1; % end point loading (N)

% analytical solution of cantilevered end point displacement
ymax = -P*L^3/3/E/I; % end displacement (m)

% mesh 
dL = L/n;

% stiffness matrix for each element
dK = [ ...
    12/dL,       6,      -12/dL,      6; ...
    6,          4*dL,    -6,         2*dL; ...
    -12/dL,      -6,     12/dL,       -6; ...
    6,          2*dL,    -6,         4*dL; ...
    ] ...
    *E*I/dL^2;

% concatenate all the elements stiffness matrices together
K = zeros(4+(n-1)*2,4+(n-1)*2);
for i = 1:n
    indx = 1+(i-1)*2;
   K(indx:indx+3,indx:indx+3) = K(indx:indx+3,indx:indx+3) + dK;
end

% cantilevered beam
F = zeros(4+(n-1)*2,1);
F(end-1) = -P;

% don't need to solve the first 2 equations due to boundary conditions
F = F(3:end);
K = K(:,3:end);
K = K(3:end,:);

% numerical solution
% x = sparse(K)\sparse(F);
% xu = x(end-1)

% our answer solution at this ID should match ymax
id = length(F)-1;


