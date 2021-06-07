%% Problem 3.4 - LP Test of Interior Point Implementation
clear
cd(fileparts(which('Problem_3_4.m')))

g = [-16.10, -8.500, -15.70, -10.02, -18.68]';
A = [   16.10,	1.0;
        8.500,	1.0;
        15.70,	1.0;
        10.02,	1.0;
        18.68,	1.0];
b = [15, 1]';

[n, m] = size(A);

l = zeros(n,1);
u = ones(n,1);

C = [eye(n) -eye(n)];
d = [l; -u];

[x, y, z, s, k] = PrimalDualInteriorLPSolver(g, A, b, C, d)
options = optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'off');
[xbar,~,~,output,lambda] = linprog(g, -C', -d, A', b, [], [], [], options);