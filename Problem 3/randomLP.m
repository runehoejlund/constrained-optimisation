function [g, A, b, C, d, x0, y0, z0] = randomLP(n,m)
    % Generates a random Linear Program
    % Inputs
    % n: number of variables
    % m: number of constraints
    %
    % Outputs
    % g, A, b, C, d of QP and solution x0, y0, z0
    
    x0 = rand(n,1);
    y0 = rand(m,1);
    z0 = rand(2*n,1);
    
    % Define g, b and d such that x0, y0, z0 solves the KKT-conditions.
    A = rand(n,m);
    b = A'*x0;
    
    l = zeros(n,1);
    u = ones(n,1);
    C = [eye(n) -eye(n)];
    d = [l; -u];
    
    g = A*y0 + C*z0;
end