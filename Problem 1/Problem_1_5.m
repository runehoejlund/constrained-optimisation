%% Problem 1.5 Equality Constrained Convex QP
clear
cd(fileparts(which('Problem_1_5.m')))
solvers = ["LUdense", "LUsparse", "LDLsparse", "LDLdense", "RangeSpace", "NullSpace"];

s = length(solvers);
maxp = 11; % Maximum problem size in power of 2.
times = nan(maxp,s);
problemSizes = nan(maxp,1);

%% Generate maxp QPs of size given in power of 2; 
for p=1:maxp
    n = 2^p;
    m = n;
    problemSizes(p) = n+m;
    [H, g, A, b, x0, lambda0] = randomQP(n,m);
    [x, lambda, times(p,:)] = testSolvers(H, g, A, b, solvers);
    
    % Assert that solution is correct.
    assert(isequal(round(x,4),round(x0,4)), " The solution should be correct");
    assert(isequal(round(lambda,4),round(lambda0,4)), " The Lagrange multipliers should be correct");
end

%% Cache the calculations for next run.
save('Problem_1_5.mat','problemSizes','times','solvers')

%% Plot figure comparing performance of different solvers
figure('Position', [0 0 500 350]);
loglog(problemSizes, times, '-o')
xlabel('Problem size, $n + m$','Interpreter','latex')
ylabel('calculation time [$s$]','Interpreter','latex')
legend(solvers, 'Location', 'northwest')
savePDF('./problem_1_5.pdf')
