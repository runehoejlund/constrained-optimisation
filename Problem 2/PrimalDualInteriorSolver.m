function [x, y, z, s, k] = PrimalDualInteriorSolver(H, g, A, b, C, d)
    [x, y, z, s] = getInitialPoint(H, g, A, b, C, d);
    
    mc = size(C,2);
    eta = 0.995;
    k = 1;
    
    % Parameters for convergence criteria
    epsilon = 10^(-6);
    mu0 = (s'*z)/mc;
    
    while true
        % Compute duality and residuals
        mc = size(C, 2);
        mu = (s'*z)/mc;
        [S, Z, rL, rA, rC, rSZ] = computeResiduals(H, g, A, b, C, d, x, y, z, s);
        
        % Check convergence
        rLCrit = norm(rL, inf) < epsilon*max(1,norm([H g A C], inf));
        rACrit = norm(rA, inf) < epsilon*max(1,norm([A' b], inf));
        rCCrit = norm(rC, inf) < epsilon*max(1,norm([eye(mc) d C'], inf));
        muCrit = mu < epsilon*10^(-2)*mu0;
        if rLCrit && rACrit && rCCrit && muCrit
            break
        end
        
        % Compute affine Newton step
        [~, ~, dzAff, dsAff, L, D, p] = NewtonStep(H,A,C,S,Z,rL,rA,rC,rSZ);
        alphaAff = computeAlpha(z, dzAff, s, dsAff);
        
        % Set centering parameter from duality measure
        zAff = z + alphaAff*dzAff;
        sAff = s + alphaAff*dsAff;
        muAff = zAff'*sAff/mc;
        sigma = (muAff/mu)^3;
        
        % Correct and center direction
        DS = diag(dsAff);
        DZ = diag(dzAff);
        e = ones(mc,1);
        rSZbar = rSZ + DS*DZ*e - sigma*mu*e;
        [dx, dy, dz, ds] = NewtonStep(H,A,C,S,Z,rL,rA,rC,rSZbar,L,D,p);
        
        % Compute alpha
        alpha = computeAlpha(z, dz, s, ds);
        
        % Update iteration
        x = x + eta*alpha*dx;
        y = y + eta*alpha*dy;
        z = z + eta*alpha*dz;
        s = s + eta*alpha*ds;
        
        k = k + 1;
    end
end

function [x, y, z, s] = getInitialPoint(H, g, A, b, C, d)
    m = size(A,2);
    [n, mc] = size(C);
    x0 = ones(n,1);
    y0 = ones(m,1);
    z0 = ones(mc,1);
    s0 = ones(mc,1);
    
    [S, Z, rL, rA, rC, rSZ] = computeResiduals(H, g, A, b, C, d, x0, y0, z0, s0);
    [dxAff, dyAff, dzAff, dsAff] = NewtonStep(H,A,C,S,Z,rL,rA,rC,rSZ);
    x = x0 + dxAff;
    y = y0 + dyAff;
    z = max(1, abs(z0 + dzAff));
    s = max(1, abs(s0 + dsAff));
end

function [S, Z, rL, rA, rC, rSZ] = computeResiduals(H, g, A, b, C, d, x, y, z, s)
    % Compute residuals
    mc = size(C,2);
    e = ones(mc,1);
    S = diag(s);
    Z = diag(z);
    rL =  H*x-A*y-C*z+g;
    rA =  -A'*x+b;
    rC =  -C'*x+s+d;
    rSZ = S*Z*e;
end

function [dx, dy, dz, ds, L, D, p] = NewtonStep(H,A,C,S,Z,rL,rA,rC,rSZ,L,D,p)
    % Compute affine direction using LDL factorisation
    % The function accepts the L D p vectors as optional arguments
    % so they can be reused if they have already been calculated
    [n, m] = size(A);
    if nargin < 10
        Hbar = H + C*(S \ Z)*C';
        Kbar = sparse([ Hbar     -A;
                    -A'   zeros(m)]);
        [L, D, p] = ldl(Kbar, 'lower', 'vector');
    end
    rLbar = rL - C*(S \ Z)*(rC - Z \ rSZ);
    dbar = - [rLbar; rA];
    sol(p,1) = L' \ ( D \ ( L \ dbar(p) ) );
    dx = sol(1:n);
    dy = sol(end-m+1:end);
    dz = S \ (-rSZ + Z*(rC  - C'*dx));
    ds = - Z \ (rSZ + S*dz);
end

function [alpha] = computeAlpha(z, dz, s, ds)
    iz = dz < 0;
    is = ds < 0;
    alpha = min([1,...
            min(- z(iz)./dz(iz)),...
            min(- s(is)./ds(is))]);
end