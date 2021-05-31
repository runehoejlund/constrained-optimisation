function [x, lambda] = EqualityQPSolverLUdense(H, g, A, b)
    [x, lambda] = EqualityQPSolver(H, g, A, b, 'LUdense');
end