function [x, lambda] = EqualityQPSolverNullSpace(H, g, A, b)
    [x, lambda] = EqualityQPSolver(H, g, A, b, 'NullSpace');
end