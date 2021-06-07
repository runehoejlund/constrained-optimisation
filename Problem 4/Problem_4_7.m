%% Problem 4 - Nonlinear Program (NLP) 
% Himmelblau's Problem - Solved with Line Search SQP + damped BFGS implementation.
clear
close all
cd(fileparts(which('Problem_4_7.m')))

% 4.4
x0s = [ 0, -4, 0.1;
        0, 0, 2.9];
xl = -5; xu = 5;

% Draw contour plot
figure(1)
fx1x2 = @(x1,x2) Himmelblau([x1; x2]);
fcontour(fx1x2,'LevelList',linspace(0,250,50))
hold on
x1 = linspace(xl,xu);
x2c1 = @(x1) (x1+2).^2;
x2c2 = @(x1) 4/10*x1;
fill([xl x1 xu],[xu min(x2c1(x1),xu) xu],[0.7 0.7 0.7],'LineWidth',1,'facealpha',0.2)
fill([xl x1 xu],[xl max(x2c2(x1),xl) xl],[0.7 0.7 0.7],'LineWidth',1,'facealpha',0.2)

for i=1:size(x0s,2)
    x0 = x0s(:,i);
    c = constraintsHimmelblau(x0);
    m = size(c,1);
    lambda0 = ones(m,1);
    lambdaForHessian.ineqnonlin = lambda0;
    B = hessianHimmelblau(x0, lambdaForHessian);
    B0 = diag(abs(eig(B)));
    [x, lambda, f, k] = SQPLineSearchSolver(@Himmelblau, @constraintsHimmelblau, x0, lambda0, B0);
    plot(x(1,:)',x(2,:)','-o','LineWidth',1);
    x, lambda, f
end
hold off

xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
savePDF('./problem_4_7.pdf')