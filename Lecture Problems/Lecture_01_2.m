%% Problem 1 from Lecture Notes: Gradient and Hessian of Multivariate Scalar Function

% 1.
clear;

fun = @(x1,x2,p) p(1)*(x2-x1.^2).^2+p(2)*(1-x1).^2;

[X1,X2] = meshgrid(-5:0.05:5,-5:0.05:20);

figure(1)
contour(X1,X2,fun(X1,X2,[100;1]),50)

%%
[f, df, d2f] = Rosenbrock([2; 3],[100;1]);

function [f, df, d2f] = Rosenbrock(x,p)
    fun = @(x) p(1)*(x(2)-x(1).^2).^2+p(2)*(1-x(1)).^2;
    dfun = @(x) [   -4*p(1)*(x(2)*x(1) - x(1)^3) - 2*p(2)*(1 - x(1));
                    2*p(1)*(x(2)-x(1)^2)];
    f = fun(x);
    df = dfun(x);
    d2f = JacobianFDforward(dfun,x);
end