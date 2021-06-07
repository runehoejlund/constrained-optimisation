%% Problem 4 - Nonlinear Program (NLP)
% Himmelblau's Problem - Contour plot
clear
close all
cd(fileparts(which('Problem_4_4.m')))

%% 4.4 Draw contour plot of Himmelblau
xl = -5; xu = 5;

figure('Position', [0 0 500 400])
fx1x2 = @(x1,x2) Himmelblau([x1; x2]);
fcontour(fx1x2,'LevelList',linspace(0,250,50))
hold on
x1 = linspace(xl,xu);

% Plot stationary points
x0s = [-3.6, -3.4, -3, -1.5, -0.45, -0.3, 0.1, 3, 3.2;
       2.5, -1.39, 0, 0.18, -0.18, 2.9, 2.88, 2, 1.3];
plot(x0s(1,:),x0s(2,:),'r.','MarkerSize',30)

% Draw constraints
x2c1 = @(x1) (x1+2).^2;
x2c2 = @(x1) 4/10*x1
fill([xl x1 xu],[xu min(x2c1(x1),xu) xu],[0.7 0.7 0.7],'LineWidth',1,'facealpha',0.2)
fill([xl x1 xu],[xl max(x2c2(x1),xl) xl],[0.7 0.7 0.7],'LineWidth',1,'facealpha',0.2)
hold off

xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
savePDF('./problem_4_4.pdf')