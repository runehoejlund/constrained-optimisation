%% Exercise 2. Problem 1: Quadratic Optimization
clear;

%% 4. Find the optimal solution
H = [6 2 1;
     2 5 2;
     1 2 4];
g = [-8; -3; -3];
A = [1 0;
     0 1;
     1 1];
b = [3; 0];
[n, m] = size(A);
K = [H  -A;
     -A' zeros(m)];
d = [-g; -b];
z0 = K \ d
x0 = z0(1:n);
lambda0 = z0(end-m+1:end);

%% 5. Sufficient conditions:
eig(H)
% All eigenvalues are positive, so H is positive definite meaning the
% problem is strictly convex and the first order conditions are therefore
% sufficient.

% The above shows that w' H w > 0 for all w (since H is positive definite). 
% We can verify that this is also true for the special case of
% the critical cone (see Important slides > Sufficient Conditions):
round(A'*x0-b,4) % all constraints are active at x0
w = null(A','r') % The null set of A' is the critical cone (here the set of feasible directions).
w'*H*w % This is positive, verifying again that the sufficient condition is met.

%% Alternative implementations for solving the problem:

% Alternatively using fmincon.
fun = @(x) 0.5*x'*H*x + g'*x;
[x0,~,~,~,lambda,~,~] = fmincon(fun,zeros(n,1),[],[],A',b);
x0
lambda.eqlin

% Alternatively using my own EqualityQPSolver.m.
[x1, lambda1] = EqualityQPSolver(H, g, A, b, "LUdense")