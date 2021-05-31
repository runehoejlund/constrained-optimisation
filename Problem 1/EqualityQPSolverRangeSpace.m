function [x, lambda] = EqualityQPSolverRangeSpace(H, g, A, b)
    [x, lambda] = EqualityQPSolver(H, g, A, b, 'RangeSpace');
end