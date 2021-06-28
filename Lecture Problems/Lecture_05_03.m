%% Lecture 5. Problem 3. Markowitz Portfolio Optimization
% 5.3

H = [   2.30 0.93 0.62 0.74 -0.23;
        0.93 1.40 0.22 0.56 0.26;
        0.62 0.22 1.80 0.78 -0.27;
        0.74 0.56 0.78 3.40 -0.56;
        -0.23 0.26 -0.27 -0.56 2.60];
g = zeros(5,1);
lb = zeros(5,1);
A = [   15.10 12.50 14.70 9.02 17.68;
        ones(1,5)];
R = 10;
b = [R 1]';
x = quadprog(H,g,[],[],A,b,lb)
%% 5.4
    
N = 100;
x = zeros(5,N);
F = zeros(1,N);
R = linspace(9.02,17.68,N);
for i = 1:N
    b = [R(i) 1]';
    [x(:,i), ~, exitflag] = quadprog(H,g,[],[],A,b,lb);
    F(:,i) = 1/2 * x(:,i)' * H * x(:,i);
    if exitflag == -2
        x(:,i) = nan(5,1);
    end
end
%%
figure(1)
subplot(2,1,1)
plot(R,F)
subplot(2,1,2)
plot(R,x)
legend('1','2','3','4','5')

%% 5.6 Risk free security
H_risk_free = [H zeros(5,1); zeros(1,5) 0]
M = size(H_risk_free,2)
g = zeros(M,1);
lb = zeros(M,1);
A = [   15.10 12.50 14.70 9.02 17.68 2;
        ones(1,M)];
    
N = 100;
x = zeros(M,N);
F = zeros(1,N);
R = linspace(9.02,17.68,N);
for i = 1:N
    b = [R(i) 1]';
    [x(:,i), ~, exitflag] = quadprog(H_risk_free,g,[],[],A,b,lb);
    F(:,i) = 1/2 * x(:,i)' * H_risk_free * x(:,i);
    if exitflag == -2
        x(:,i) = nan(M,1);
    end
end

%%
figure(2)
subplot(1,3,1)
plot(R,F)
xline(15)
subplot(1,3,2)
plot(R,x)
legend('1','2','3','4','5','6')
xline(15)
subplot(1,3,3)
[~, minidx] = min(abs(R-15));
pie3(x(:,minidx),[0,1,1,2,3,4])