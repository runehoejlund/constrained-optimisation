%% Problem 4 - Nonlinear Program (NLP)
% Himmelblau's Problem - Solved with fmincon in Matlab
clear
close all
cd(fileparts(which('Problem_4_5.m')))

%% 4.5 Solve Himmelblaus probelm with fmincon and by providing gradients
x0s = [ 0, -4, 0.1;
        0, 0, 2.9];
options = optimoptions('fmincon',...
        'Algorithm','interior-point',...
        "SpecifyConstraintGradient",true,...
        "SpecifyObjectiveGradient",true,...
        'Display','none');

xl = -5; xu = 5;
xs = nan(size(x0s));
for i=1:length(x0s)
    x0 = x0s(:,i);
    xs(:,i) = fmincon(@Himmelblau,x0,[],[],[],[],[],[],@fminconstraintsHimmelblau,options,[]);
end

%% Draw contour plot of solutions and starting points
figure('Position', [0 0 500 400])
fx1x2 = @(x1,x2) Himmelblau([x1; x2]);
fcontour(fx1x2,'LevelList',linspace(0,250,50))
hold on
x1 = linspace(xl,xu);

% Draw constraints
x2c1 = @(x1) (x1+2).^2;
x2c2 = @(x1) 4/10*x1;
fill([xl x1 xu],[xu min(x2c1(x1),xu) xu],[0.7 0.7 0.7],'LineWidth',1,'facealpha',0.2)
fill([xl x1 xu],[xl max(x2c2(x1),xl) xl],[0.7 0.7 0.7],'LineWidth',1,'facealpha',0.2)

% Plot starting points
plot(x0s(1,:),x0s(2,:),'bx','MarkerSize',10,'LineWidth',2)
plot(xs(1,:),xs(2,:),'r.','MarkerSize',30)

hold off

xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
savePDF('./problem_4_5.pdf')