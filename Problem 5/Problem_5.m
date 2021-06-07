%% Problem 5. Markowitz Portfolio Optimization
% 5.3
clear
close all
cd(fileparts(which('Problem_5.m')));

options = optimoptions('quadprog',...
    'Algorithm','interior-point-convex',...
    'Display','none');

H = [   2.30 0.93 0.62 0.74 -0.23;
        0.93 1.40 0.22 0.56 0.26;
        0.62 0.22 1.80 0.78 -0.27;
        0.74 0.56 0.78 3.40 -0.56;
        -0.23 0.26 -0.27 -0.56 2.60];
n = size(H,1);
e = ones(n,1);
g = zeros(5,1);
lb = zeros(5,1);
A = [   15.10 12.50 14.70 9.02 17.68;
        e']';
Rs = 10;
b = [Rs 1]';
x = quadprog(H,g,[],[],A',b,lb)
0.5*x'*H*x

%% 5.4
    
N = 100;
xs = zeros(5,N);
Fs = zeros(1,N);
Rs = linspace(9.02,17.68,N);
for i = 1:N
    b = [Rs(i) 1]';
    [xs(:,i), ~, ~] = quadprog(H,g,[],[],A',b,lb,[],[],options);
    Fs(:,i) = 0.5 * xs(:,i)' * H * xs(:,i);
end

figure(1)
subplot(2,1,1)
plot(Rs,Fs)
ylabel('Risk $F$','Interpreter','latex')
subplot(2,1,2)
plot(Rs,xs)
legend('Asset 1','Asset 2','Asset 3','Asset 4','Asset 5')
xlabel('Return, $R$','Interpreter','latex')
ylabel('Distribution $x_i$','Interpreter','latex')
savePDF('./problem_5_4.pdf')

%% 5.5 Risk free security
H_risk_free = [H zeros(5,1); zeros(1,5) 0]
M = size(H_risk_free,2);
e = ones(M,1);
g = zeros(M,1);
lb = zeros(M,1);
A = [   15.10 12.50 14.70 9.02 17.68 0;
        e']';
    
N = 100;
xs = zeros(M,N);
Fs = zeros(1,N);
Rs = linspace(9.02,17.68,N);
for i = 1:N
    b = [Rs(i) 1]';
    [xs(:,i), ~, exitflag] = quadprog(H_risk_free,g,[],[],A',b,lb,[],[],options);
    Fs(:,i) = 1/2 * xs(:,i)' * H_risk_free * xs(:,i);
end

%%
figure(2)
subplot(2,1,1)
plot(Rs,Fs)
ylabel('Risk $F$','Interpreter','latex')
xline(14,'-','$R = 14$','Interpreter','latex')
subplot(2,1,2)
plot(Rs,xs)
xline(14,'-','$R = 14$','Interpreter','latex')
legend('Asset 1','Asset 2','Asset 3','Asset 4','Asset 5', 'Risk Free')
xlabel('Return, $R$','Interpreter','latex')
ylabel('Distribution $x_i$','Interpreter','latex')
savePDF('./problem_5_risk_free.pdf')

%% 
Rs = 14;
b = [Rs 1]';
x = quadprog(H_risk_free,g,[],[],A',b,lb)
0.5*x'*H_risk_free*x