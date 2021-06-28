%% Problem 1 from Lecture Notes: Gradient and Hessian of Multivariate Scalar Function

% 1.
clear;

fun = @(x1,x2) x1.^2 - 2*x1 + 3*x1.*x2 + 4*x2.^3;

[X1,X2] = meshgrid(-5:0.01:5,-5:0.01:5);

figure(1)
contour(X1,X2,fun(X1,X2),50)

%%
[f, df, d2f] = FunEx1([2; 3])
[f, df, d2f] = FunEx1Alternative([2; 3])

function [f, df, d2f] = FunEx1(x)
    fun = @(x) x(1).^2 - 2*x(1) + 3*x(1).*x(2) + 4*x(2).^3;
    [f, df, d2f] = Differentiate(fun,x);
end

%% Alternative
% if we can use the analytical expression for the gradient:

function [f, df, d2f] = FunEx1Alternative(x)
    fun = @(x) x(1).^2 - 2*x(1) + 3*x(1).*x(2) + 4*x(2).^3;
    dfun = @(x) [2*x(1) - 2 + 3 *x(2); 3*x(1) + 12*x(2)^2];
    f = fun(x);
    df = dfun(x);
    d2f = JacobianFDforward(dfun,x);
end