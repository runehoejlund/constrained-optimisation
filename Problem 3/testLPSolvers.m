function [x, y, z, s, ks, times, solvers] = testLPSolvers(g, A, b, C, d)
    % Test the different solvers
    solvers = ["My Implementation", "linprog"];
    s = length(solvers);
    [n, m] = size(A);
    mc = size(C,2);
    times = nan(1,s);
    ks = nan(1,s);
    xs = nan(n,s);
    ys = nan(m,s);
    zs = nan(mc,s);
    ss = nan(mc,s);
    for i = 1:s
        solver = solvers(i);
        if solver == "linprog"
            options = optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'off');
            times(:,i) = timeit(@() linprog(g, -C', -d, A', b, [], [], [], options));
            [xs(:,i),~,~,output,lambda] = linprog(g, -C', -d, A', b, [], [], [], options);
            ys(:,i) = lambda.eqlin;
            zs(:,i) = lambda.ineqlin;
            ss(:,i) = C'*xs(:,i)-d;
            ks(:,i) = output.iterations;
        else
            times(:,i) = timeit(@() PrimalDualInteriorLPSolver(g, A, b, C, d));
            [xs(:,i), ys(:,i), zs(:,i), ss(:,i), ks(:,i)] = PrimalDualInteriorLPSolver(g, A, b, C, d);
        end
    end
    % Display the solution and calculation times
    x = xs(:,1);
    y = ys(:,1);
    z = zs(:,1);
    s = ss(:,1);
    
    % Tests: Assert that the different solvers agree.
    if rank(round(xs,3)) > 1 || rank(round(ss,3)) > 1
        disp(xs)
        disp(ss)
        warning("solvers disagree on solution")
    end
    primalFeasibility = round(sum(A'*xs-b),3) == 0;
    if any(~primalFeasibility)
        warning(strcat(solvers(~primalFeasibility), " gives inaccurate x for problem size n + m + mc= ", string(n + m + mc)))
    end
    
end