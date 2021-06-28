%% Problem 1: Gradient and Hessian of Multivariate Scalar Function
clc;
clear all;
x1 = -10:0.01:10;
x2 = -10:0.01:10;
[X1,X2] = meshgrid(x1,x2);
f = @(X1, X2) (X1.^2 - 2*X1 + 3.*X1.*X2 + 4*X2.^3);

figure(1)
contour(X1,X2,f(X1,X2))

%% 4.
[f, df, d2f] = FunEx1([2, 3])

function [f, df, d2f] = FunEx1(x)
    fun = @(x1,x2) (x1.^2 - 2*x1 + 3.*x1.*x2 + 4*x2.^3);
    f = fun(x(1),x(2));
    e = sqrt(eps);
    dfun = @(x1,x2) ([fun(x1 + e,x2) - fun(x1,x2); fun(x1, x2 + e) - fun(x1,x2)]./e); 
    df = dfun(x(1),x(2));
    d2f = nan(1);
    %% BAD ATTEMPT
    %     ddfun = @(x1,x2) [dfun(x1 + e,x2) - dfun(x1,x2); dfun(x1, x2 + e) - dfun(x1,x2)]./e;
    %     d2f = ddfun(x(1),x(2));
end