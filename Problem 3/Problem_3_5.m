%% Problem 3.5 - LP Test of Interior Point Implementation
clear
cd(fileparts(which('Problem_3_5.m')))
solvers = ["My Implementation", "linprog"];

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
    [g, A, b, C, d] = randomLP(n,m);
    
    [x, y, z, s, ks, times(p,:)] = testLPSolvers(g, A, b, C, d);
end

%% Cache the calculations for next run.
save('Problem_3_5.mat','problemSizes','times','solvers')

%% Plot figure comparing performance of different solvers
figure('Position', [0 0 500 350]);
loglog(problemSizes, times, '-o')
xlabel('Problem size, $n + m$','Interpreter','latex')
ylabel('calculation time [$s$]','Interpreter','latex')
legend(solvers, 'Location', 'northwest')
savePDF('./problem_3_5.pdf')
