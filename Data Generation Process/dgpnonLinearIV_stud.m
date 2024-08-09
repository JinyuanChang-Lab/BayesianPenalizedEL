function [y, X, Z] = dgpnonLinearIV_stud(b, rho, useful)

%rho: magnitude of endogeneity
global n m 

%% instruments and disturbances

Z = random('T', 3, n, m );
E = random('Normal', 0, 1, n, 3);
Sigma = [ 1, rho, rho; rho, 1, 0; rho, 0, 1 ];
E = E * Sigma;
e = E(:,1);
v1 = E(:,2);
v2 = E(:,3);

%% structural equation
uu = useful/2;

x1 = 1/sqrt( uu ) * Z(:, 1: uu )  * ones(uu, 1) + 0.5 * v1;
x2 = 1/sqrt( uu ) * Z(:, (uu+1):useful )  * ones(uu, 1) + 0.5 * v2;

X = [x1, x2];

y = sin(X * b) + 0.5 * e;

end
