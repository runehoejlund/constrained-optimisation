%% Problem 2.4 - QP Test of Interior Point Implementation
clear
cd(fileparts(which('Problem_2_5.m')))
solvers = ["My Implementation", "quadprog"];

s = length(solvers);
maxp = 11; % Maximum problem size in power of 2.
times = nan(maxp,s);
problemSizes = nan(maxp,1);
ks = nan(maxp,s);

%% Generate maxp QPs of sizes given in powers of 2; 
for p=1:maxp
    n = 2^p;
    m = n/2;
    problemSizes(p) = n+m;
    [H, g, A, b, C, d, x0, y0, z0] = randomIQP(n,m);
    [x, y, z, s, ks, times(p,:)] = testQPSolvers(H, g, A, b, C, d, solvers);
end

%% Cache the calculations for next run.
save('Problem_2_5.mat','problemSizes','times','solvers')

%% Plot figure comparing performance of different solvers
figure('Position', [0 0 500 350]);
loglog(problemSizes, times, '-o')
xlabel('Problem size, $n + m$','Interpreter','latex')
ylabel('calculation time [$s$]','Interpreter','latex')
legend(solvers, 'Location', 'northwest')
savePDF('./problem_2_5.pdf')
