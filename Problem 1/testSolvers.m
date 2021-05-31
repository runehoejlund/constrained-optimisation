function [x, lambda, times] = testSolvers(H, g, A, b, solvers)
    % Test the different solvers
    s = length(solvers);
    [n, m] = size(A);
    times = nan(1,s);
    xs = nan(n,s);
    lambdas = nan(m,s);
    for i = 1:s
        solver = solvers(i);
        times(:,i) = timeit(@() EqualityQPSolver(H, g, A, b, solver));
        [xs(:,i), lambdas(:,i)] = EqualityQPSolver(H, g, A, b, solvers(i));
    end
    % Display the solution and calculation times
    x = xs(:,1);
    lambda = lambdas(:,1);
    
    % Tests: Assert that the different solvers agree.
    if rank(xs,1e-3) ~= 1 || rank(lambdas,1e-3) ~= 1
        warning("solvers disagree on solution")
    end
    primalFeasibility = round(sum(A'*xs-b),3) == 0;
    if any(~primalFeasibility)
        warning(strcat(solvers(~primalFeasibility), " gives inaccurate x for problem size n + m = ", string(n+m)))
    end
    stationarity = round(sum(A*lambdas-H*xs-g),3) == 0;
    if any(~stationarity)
        warning(strcat(solvers(~stationarity), " gives inaccurate lambda for problem size n + m = ", string(n+m)))
    end
end