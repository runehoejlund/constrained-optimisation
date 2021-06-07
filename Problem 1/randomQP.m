function [H, g, A, b, x0, lambda0] = randomQP(n,m)
    % Generates a random Convex Quadratic Problem
    % Inputs
    % n: number of variables
    % m: number of constraints
    %
    % Outputs
    % H, g, A, b, C, d of QP and solutions x0 and lambda0 and z0 (multiplier).
    
    x0 = rand(n,1);
    lambda0 = rand(m,1);
    
    % Define squareroot of H so that H is a positive operator.
    sqrtH = rand(n);
    H = sqrtH'*sqrtH;
    A = rand(n,m);
    b = A'*x0;
    
    % Define g and b such that x0 and lambda0 (and z0) solves the KKT-conditions.
    g = A*lambda0 - H*x0;
end