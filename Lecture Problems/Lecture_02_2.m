%% Exercise 2. Problem 2: Linear Optimization
clear;
close all;

% 2.1
g = [1; -2];
A = [1  0   1   1   -5
     0  1   -1  -5  1];
b = [0; 0;  -2; -20;    -15];
f = @(x1, x2) g(1)*x1+g(2)*x2;

x1 = -1:0.01:5;
x2 = -1:0.01:5;

figure(1)
[X1,X2] = meshgrid(x1,x2);
contour(X1,X2,f(X1,X2),-15:1:15,'ShowText','on');
plotConstraints(A,b,x1,x2)

%% 2.5
% {3, 4} is the active set
% so the equations to solve is
% [a3 a4] [lambda3; lambda4] = g and [a3 a4]'*x = [b3; b4]
Aactive = [A(:,3) A(:,4)]
bactive = [b(3); b(4)]
x = Aactive' \ bactive
lambda34 = Aactive \ g
lambda = [0; 0; lambda34; 0]

% Check
assert(all(g-A*lambda == 0), "stationarity")
assert(all(A'*x-b >= 0), "feasibility")