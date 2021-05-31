function [x, lambda] = EqualityQPSolverLDLsparse(H, g, A, b)
    [x, lambda] = EqualityQPSolver(H, g, A, b, 'LDLsparse');
end