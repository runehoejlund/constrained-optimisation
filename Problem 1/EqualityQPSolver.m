function [x, lambda] = EqualityQPSolver(H, g, A, b, solver)
    K = sparse([ H     -A;
          -A'   zeros(size(A,2))]);
    d = - [g; b];
    switch solver
        case 'LUdense'
            [L, U, p] = lu(K, 'vector');
            z = U \ ( L \ d(p) );
        case 'LUsparse'
            [L, U, p, q] = lu(K, 'vector');
            z(q,1) = U \ ( L \ d(p) );
        case 'LDLdense'
            [L, D, p] = ldl(K, 'lower', 'vector');
            z(p,1) = L' \ ( D \ ( L \ d(p) ) );
        case 'LDLsparse'
            [L, D, P, S] = ldl(K, 'lower');
            z = S * P * (L' \ ( D \ ( L \  (P' * S * d) )));
        case 'RangeSpace'
            invH = inv(H);
            lambda = (A' * invH * A) \ (b + A' * invH * g);
            x = invH * (A * lambda - g);
            z = [x; lambda];
        case 'NullSpace'
            Z = null(A');
            Y = null(Z');
            xY = (A'*Y) \ b;
            xZ = - (Z'*H*Z) \ (Z' * (H*Y*xY + g));
            x = Y*xY + Z*xZ;
            lambda = (A'*Y)' \ (Y' * (H*x + g));
            z = [x; lambda];
    end
    [n, m] = size(A);
    x = z(1:n); % n is dimension of x.
    lambda = z(end-m+1:end); % m is number of constraints.
end