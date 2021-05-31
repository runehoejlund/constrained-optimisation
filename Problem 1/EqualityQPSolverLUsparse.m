function [x, lambda] = EqualityQPSolverLUsparse(H, g, A, b)
    [x, lambda] = EqualityQPSolver(H, g, A, b, 'LUsparse');
end