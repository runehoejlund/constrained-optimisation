%% Problem 1.4 - Equality Constrained Convex QP
clear
cd(fileparts(which('Problem_1_4.m')))

solvers = ["LUdense", "LUsparse", "LDLsparse", "LDLdense", "RangeSpace", "NullSpace"];

H = [   5.000,	1.86,	1.240,	1.480,	-0.46;
        1.860,	3.00,	0.440,	1.120,	0.520;
        1.240,	0.44,	3.800,	1.560,	-0.54;
        1.480,	1.12,	1.560,	7.200,	-1.12;
        -0.46,	0.52,	-0.54,	-1.12,	7.800];
g = [-16.10, -8.500, -15.70, -10.02, -18.68]';
A = [   16.10,	1.0;
        8.500,	1.0;
        15.70,	1.0;
        10.02,	1.0;
        18.68,	1.0];
b = [15, 1]';

[x, lambda, times] = testSolvers(H, g, A, b, solvers)

% Test that the solvers work
assert(all(round(A'*x-b,3) == 0), " The KKT-conditions shoud be satisfied")
assert(all(round(H*x-A*lambda + g,3) == 0), " The KKT-conditions should be satisfied")

%%
% Calculate and save solutions for different values of b1.
b1s = linspace(8.5,18.68,10);
xs = nan(size(A,1),length(b1s));
lambdas = nan(size(A,2),length(b1s));
for i=1:length(b1s)
    b(1) = b1s(i);
    [xs(:,i), lambdas(:,i), times] = testSolvers(H, g, A, b, solvers);
    assert(all(round(A'*xs(:,i)-b,3) == 0), " The KKT-conditions shoud be satisfied")
    assert(all(round(H*xs(:,i)-A*lambdas(:,i) + g,3) == 0), " The KKT-conditions should be satisfied")
end

%% Plot figure comparing solution for different values of b1
figure('Position', [0 0 500 250]);
plot(b1s,xs,'o')
legend(strcat("$x_",string(1:5),"$"),'Interpreter','latex')
xlabel('$b_1$','Interpreter','latex')
ylabel('Coordinates of minimizer','Interpreter','latex')
savePDF('./problem_1_4.pdf')