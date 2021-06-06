function [x, y, z, s, k] = PrimalDualInteriorSolver(H, g, A, b, C, d, x0, y0, z0, s0)
    assert(all(z0 > 0), "Lagrange multipliers must be positive: z0 > 0")
    assert(all(s0 > 0), "Inequality constraints must be positive: s0 > 0")
    x = x0; y = y0; z = z0; s = s0;
    m = size(A,2);
    [n, mc] = size(C);
    eta = 0.995;
    e = ones(mc,1);
    epsilon = 10^(-3);
    mu0 = (s'*z)/mc;
    
    k = 0;
    
    while true
        % Duality and residuals
        mu = (s'*z)/mc;
        S = diag(s);
        Z = diag(z);
        rL =  H*x-A*y-C*z+g;
        rA =  -A'*x+b;
        rC =  -C'*x+s+d;
        rSZ = S*Z*e;
        
        % Check convergence
        rLCrit = norm(rL, inf) < epsilon*max(1,norm([H g A C], inf));
        rACrit = norm(rA, inf) < epsilon*max(1,norm([A' b], inf));
        rCCrit = norm(rC, inf) < epsilon*max(1,norm([eye(mc) d C'], inf));
        muCrit = mu < epsilon*10^(-2)*mu0;
        if rLCrit && rACrit && rCCrit && muCrit
            break
        end
        
        % Compute affine direction using LDL factorisation
        Hbar = H + C*(S \ Z)*C';
        rLAff = rL - C*(S \ Z)*(rC - Z \ rSZ);
        Kbar = sparse([ Hbar     -A;
                        -A'   zeros(m)]);
        dAff = - [rLAff; rA];
        [L, D, p] = ldl(Kbar, 'lower', 'vector');
        solAff(p,1) = L' \ ( D \ ( L \ dAff(p) ) );
        dxAff = solAff(1:n);
        % [dxAff, dyAff] = EqualityQPSolver(Hbar, rLbar, A, rA, solver);
        
        dzAff = S \ (-rSZ + Z*(rC  - C'*dxAff));
        dsAff = - Z \ (rSZ + S*dzAff);
        
        % Compute affine alpha
        izAff = dzAff < 0;
        isAff = dsAff < 0;
        alphaAff = max(1, min(  min(- z(izAff)./dzAff(izAff)),...
                                min(- s(isAff)./dsAff(isAff)) ));
        
        % Duality gap and centering parameter
        zAff = z + alphaAff*dzAff;
        sAff = s + alphaAff*dsAff;
        muAff = zAff'*sAff/mc;
        sigma = (muAff/mu)^3;
        
        % Affine Centering-Correction Direction
        DS = diag(dsAff);
        DZ = diag(dzAff);
        rSZbar = rSZ + DS*DZ*e - sigma*mu*e;
        rlbar = rL - C*(S \ Z)*(rC - Z \ rSZbar);
        dbar = - [rlbar; rA];
        sol(p,1) = L' \ ( D \ ( L \ dbar(p) ) );
        dx = sol(1:n);
        dy = sol(end-m+1:end);
        dz = S \ (-rSZbar + Z*(rC  - C'*dx));
        ds = - Z \ (rSZbar + S*dz);
        
        % Compute alpha
        iz = dz < 0;
        is = ds < 0;
        alpha = eta*max(1, min(  min(- z(iz)./dz(iz)),...
                                min(- s(is)./ds(is)) ));
        
        % Update iteration
        x = x + alpha*dx;
        y = y + alpha*dy;
        z = z + alpha*dz;
        s = s + alpha*ds;
        
        k = k + 1;
    end
end
% end