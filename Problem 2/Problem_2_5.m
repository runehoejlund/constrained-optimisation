%% Problem 1.4 - Equality Constrained Convex QP
clear
cd(fileparts(which('Problem_2_5.m')))

H = [   5.000,	1.86,	1.240,	1.480,	-0.46;
        1.860,	3.00,	0.440,	1.120,	0.520;
        1.240,	0.44,	3.800,	1.560,	-0.54;
        1.480,	1.12,	1.560,	7.200,	-1.12;
        -0.46,	0.52,	-0.54,	-1.12,	7.800];
g = [-16.10, -8.500, -15.70, -10.02, -18.68]';
A = [   16.10,	1.0;
        8.500,	1.0;
        15.70,	1.0;
        10.02,	1.0;
        18.68,	1.0];
b = [15, 1]';

[n, m] = size(A)

l = zeros(n,1);
u = ones(n,1);

C = [eye(n) -eye(n)];
d = [l; -u];


%%
[x, y, z, s, k] = PrimalDualInteriorSolver(H, g, A, b, C, d)
options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex', 'Display', 'off');
[xquad,~,~,~,lambda0] = quadprog(H,g,[],[],A',b,l,u,[],options)
            
quadprog(H, g, -C', -d, A', b, [], [], [], options)

%%
[x, y, z, s, ks, times, solvers] = testQPSolvers(H, g, A, b, C, d);

% Calculate and save solutions for different values of b1.
b1s = linspace(8.5,18.68,10);
s = length(solvers);
n = size(A,1);
xs = nan(length(b1s),n);
times = nan(length(b1s),s);
ks = nan(length(b1s),s);
for i=1:length(b1s)
    b(1) = b1s(i);
    [xs(i,:), ~, ~, ~, ks(i,:), times(i,:)] = testQPSolvers(H, g, A, b, C, d);
end

%%
plot(b1s, times)
