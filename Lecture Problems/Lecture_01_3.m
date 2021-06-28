%% Problem 1 from Lecture Notes: Gradient and Hessian of Multivariate Scalar Function

% 1.
clear;

%%
[c, dc, d2c] = FunEx2([1.7; 2.1])
[c, J] = FunJacEx2([1.7; 2.1])

function [c, dc, d2c] = FunEx2(x)
    fun = @(x) [exp(x(1))-x(2); x(1)^2 - 2*x(2)];
    [c, dc, d2c] = Differentiate(fun,x);
end

function [c, J] = FunJacEx2(x)
    [c, dc] =  FunEx2(x);
    J = dc';
end