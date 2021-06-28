%% Exercise 3. Problem 3: Nonlinear Optimization
clear;
close all;

% 3.1
f = @(x1, x2) x1.^2 + x2.^2 + 3*x2;
c = @(x1, x2) x1.^2 + (x2+1).^2 -1;
x1 = -2:0.01:2;
x2 = -2:0.01:2;

figure(1)
[X1,X2] = meshgrid(x1,x2);
contour(X1,X2,f(X1,X2),'ShowText','on');
hold on
fimplicit(c,'--','LineWidth',2)

%% 3.2
figure(2)
plot(x2,sqrt(1-(x2+1).^2))
hold on
plot(x2,-sqrt(1-(x2+1).^2))

% The function is symmetric around x1 = 0
figure(3)
f1D = @(x2) f(sqrt(1-(x2+1).^2),x2);
plot(-2:0.01:0,f1D(-2:0.01:0))

%% Show that 3.3. (Double check)
eig([2  0; 
     0   2])

%% Show that 3.4
% Evident from ∂f1D/∂x2 (0) = 1.
% So we can decrease f by decreasing x2

% Also evident since grad L (0,lambda) = diag(2 [1-lambda])
% so with lambda = 3/2 the gradient of the Lagrangian
% is negative definite. Thus the 2nd order condition
% is NOT satisfied.


%% 3.5 Find the minimizer...
% x = (0; -2) is the minimizer with lambda = 1/2 > 0
% With these values:
% grad L = 0
% c(x) = 0
% and hessian of L = [1 0; 0 1] which is positive definite.
% These are the sufficient conditions