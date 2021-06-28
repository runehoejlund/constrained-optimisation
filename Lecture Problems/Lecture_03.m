%% 1. Quadratic Optimization Sensitivity and Duality
clear;
close all
H = [ 6 2 1; 2 5 2; 1 2 4];
g = [-8; -3; -3];
A = [1 0; 0 1; 1 1];
b = [3; 0];

%% 1.3 + 1.4
[n, m] = size(A);
K = [ H  -A ;
      -A' zeros(m)];

[x, lambda] = EqualityQPSolver(H, g, A, b,'LUdense')

%% 1.6
% See Notes "Numerical Methods Constrained Optimization Book.pdf" last section on sensitivity
% When using the sensitivity equations for a QP
% we simply end up with -inv(K) as the sensitivity (dz/dp)
dzdp = sensitivityQP(H,A)

%% Test it by changing b a little
dzdb = dzdp(n+1:end,:);

epsilon = 0.1;
btilde = b + epsilon*[1; 0];
[2 ; -1; 1; 3; -2] + dzdb'*(btilde - b)

% It matches perfectly with the exact result from EqualityQPSolver:
[xtilde, lambdatilde] = EqualityQPSolver(H, g, A, btilde,'LUdense')

%% Also test it by changing g a little
dzdg = dzdp(1:n,:);
gtilde = g + epsilon*[1; 1; 1];
[2 ; -1; 1; 3; -2] + dzdg'*(gtilde - g)

% It matches perfectly with the exact result from EqualityQPSolver:
[xtilde, lambdatilde] = EqualityQPSolver(H, gtilde, A, b,'LUdense')

%% Function for calculating sensitivity
function [dzdp] = sensitivityQP(H, A)
    [~, m] = size(A);
    K = [H  -A ;
        -A' zeros(m)];
    dzdp = -inv(K);
end