%% Execise Lecture 4 - 1. Use fmincon in Matlab’s Optimization Toolbox
clear
close all

% 1.a.
x0 = [0;0];
A = []; b = []; Aeq = []; beq = [];
lb = []; ub = [];
xl = [-5;-5]; xu = [5;5];
options = optimoptions('fmincon',...
        'Algorithm','interior-point',...
        'Display','none');
solve = @() fmincon(@objfunHimmelblau,x0,A,b,Aeq,beq,lb,ub,@con,options);
solve()
disp('timeit(solve)')
timeit(solve)

%% Draw contour plot
figure(1)
fx1x2 = @(x1,x2) objfunHimmelblau([x1; x2]);
fcontour(fx1x2,'LevelList',linspace(0,250,50))

%% 1.b. with gradients
p = [];
optionsGrad = optimoptions('fmincon',...
        'Algorithm','interior-point',...
        "SpecifyConstraintGradient",true,...
        "SpecifyObjectiveGradient",true,...
        'Display','none');
gradSolve = @() fmincon(@objfunGradHimmelblau,x0,A,b,Aeq,beq,lb,ub,@conGrad,optionsGrad,p)
gradSolve()
disp('timeit(gradSolve)')
timeit(gradSolve)

%% 1.b. with Hessian
optionsHessian = optimoptions('fmincon',...
        'Algorithm','interior-point',...
        "SpecifyConstraintGradient",true,...
        "SpecifyObjectiveGradient",true,...
        "HessianFcn",@hessian,...
        'Display','none');
hessSolve = @() fmincon(@objfunGradHimmelblau,x0,A,b,Aeq,beq,lb,ub,@conGrad,optionsHessian);
hessSolve()
disp('timeit(hessSolve)')
timeit(hessSolve)

%%
disp('solve | gradSolve | hessSolve')
[timeit(solve), timeit(gradSolve), timeit(hessSolve)]

%%
function [f] = objfunHimmelblau(x,p)
    f = (x(1).^2 + x(2) - 11).^2 + (x(1) + x(2).^2 - 7).^2;
end

function [f, gradf] = objfunGradHimmelblau(x,p)
    f = objfunHimmelblau(x);
    % f = (x(1).^2 + x(2) - 11).^2 + (x(1) + x(2).^2 - 7).^2;
    gradf = [4*(x(1)^2+x(2)-11)*x(1) + 2*(x(1)+x(2)^2-7);
            2*(x(1)^2+x(2)-11) + 4*(x(1)+x(2)^2-7)*x(2)];
end

function [c, ceq] = con(x)
    c = [-(x(1)+2).^2 + x(2);
        4*x(1)-10*x(2)];
    ceq = [];
end

function [c, ceq, gradc, gradceq] = conGrad(x, p)
    [c, ceq] = con(x);
    gradc = [-2*(x(1)+2)      4;
            1              -10];
    gradceq = [];
end

function h = hessian(x,lambda)
    % Hessian of f 
    hessf = [8*x(1)^2 + 4*(x(1)^2+x(2)-11) + 2, 4*x(1)+4;
         4+4*x(2),                          2+8*x(2)^2+4*(x(1)+x(2)^2-7)];
    % Hessian of c
    hessc1 = [2,    0;
              0,    0];
    hessc2 = zeros(2);
    h = hessf + lambda.ineqnonlin(1)*hessc1  + lambda.ineqnonlin(2)*hessc2;
end
