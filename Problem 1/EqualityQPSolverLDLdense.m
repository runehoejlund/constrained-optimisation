function [x, lambda] = EqualityQPSolverLDLdense(H, g, A, b)
    [x, lambda] = EqualityQPSolver(H, g, A, b, 'LDLdense');
end