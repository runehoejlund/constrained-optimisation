function [H, g, A, b, C, d, x0, y0, z0] = randomIQP(n,m)
    % Generates a random Convex Quadratic Program with inequalities
    % Inputs
    % n: number of variables
    % m: number of constraints
    %
    % Outputs
    % H, g, A, b, C, d of QP and solutions x0 and y0 and z0.
    
    x0 = rand(n,1);
    y0 = rand(m,1);
    z0 = rand(m,1);
    
    % Define squareroot of H so that H is a positive operator.
    sqrtH = rand(n);
    H = sqrtH'*sqrtH;
    
    % Define g and b and d such that x0 and y0 and z0 solves the KKT-conditions.
    A = rand(n,m);
    b = A'*x0;
    
    C = rand(n,m);
    d = C'*x0;
    
    g = A*y0 - C*z0 - H*x0;
end