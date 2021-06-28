%% Problem 1: Inequality Constrained Quadratic Program- ming
clc;
clear all;
close all;
H = [2 0; 0 2];
g = [-2; -5];
A = [1  -1  -1  1   0
    -2  -2  2   0   1];
b = [2; 6;  2;  0;  0];
x1 = -1:0.01:4;
x2 = -1:0.01:4;
[X1,X2] = meshgrid(x1,x2);

%% 1.1
f = @(X1, X2) (X1 - 1).^2 + (X2 - 2.5).^2;

figure(1)
contourf(X1,X2,f(X1,X2),100)
hold on
plotregion(A',-b)
axis equal;

%% 1.2
x = [0; 0];
Wk = [];
I = 1:length(b);
for k = 1:100
    [p, lambda] = EqualityQPSolver(H, H*x(:,k)+g, A(:,Wk), zeros(size(Wk)))
    if norm(round(p,10)) == 0
        if sum(lambda >= 0) == length(lambda)
            x(:,k)
            break;
        else
            [~, idx] = min(lambda);
            x(:,k+1) = x(:,k);
            Wk(idx) = []
        end
    else
        Ik = setdiff(I, Wk);
        Ik = Ik(A(:,Ik)'*p < 0);
        ai = A(:,Ik);
        bi = b(Ik);
        ci = ai'*x(:,k) + bi;
        [alpha, idx] = min(- ci./(ai'*p));
        l = Ik(idx);
        if alpha < 1
            x(:,k+1) = x(:,k) + alpha * p
            Wk = [Wk; l];
        else
            x(:,k+1) = x(:,k) + p
            Wk = Wk;
        end
    end
end

plot(x(1,:),x(2,:),'-or')