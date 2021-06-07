function [x, lambda, f, k] = SQPLineSearchSolver(fun, constraints, x0, lambda0, B0)
    B = B0;
    k = 1;
    x(:,k) = x0;
    lambda(:,k) = lambda0;
    mu = 2*abs(lambda0);
    [f(k), df(:,k)] = fun(x(:,k));
    [c(:, k), dc] = constraints(x(:,k));
    
    % Parameters for convergence criteria
    epsilon = 10^(-6);
    
    % Parameters for Line Search
    eta = 0.1;
    
    while true
        % Check convergence
        if norm(df(:,k) - dc*lambda(:,k), inf) < epsilon
            break
        end
        
        % Compute x_(k+1) and lambda_(k+1)
        options = optimoptions('quadprog',...
            'Algorithm','interior-point-convex',...
            'Display','none');
        [p,~,~,~,lambdas] = quadprog(B, df(:,k), -dc', c(:,k), [],[], [], [], [], options);
        lambdahat = lambdas.ineqlin;
        % Line Search
        plambda = lambdahat - lambda(:,k);
        mu = max(abs(lambdahat), 0.5*(mu + abs(lambdahat)));
        alpha = 1;
        phi0 = f(k) + mu'*abs(min(0,c(:,k)));
        dphi0 = df(:,k)'*p - mu'*abs(min(0,c(:,k)));
        while true
            x1 = x(:,k) + alpha*p;
            f1 = fun(x1);
            c1 = constraints(x1);
            phi = f1 + mu'*abs(min(0,c1));
            if phi <= phi0 + eta*alpha*dphi0
                break
            else
                a = (phi - (phi0 + alpha*dphi0))/(alpha^2);
                alphamin = -dphi0/(2*a);
                alpha = min(0.9*alpha, max(alphamin,0.1*alpha));
            end
        end
        x(:,k+1) = x(:,k) + alpha*p;
        lambda(:, k+1) = lambda(:,k) + alpha*plambda;
        
        % Compute gradient of Lagrangian in x_k
        dLk = df(:,k) - dc*lambda(:,k+1);
        
        % Evaluate objective and constraints in x_(k+1)
        [f(k+1), df(:,k+1)] = fun(x(:,k+1));
        [c(:, k+1), dc] = constraints(x(:,k+1));
        
        % Compute gradient of Lagrangian in x_(k+1)
        dLkp1 = df(:,k+1) - dc*lambda(:,k+1);
        q = dLkp1 - dLk;
        
        % Damped BFGS Updating
        B = BFGS(B,alpha*p,q);
        k = k + 1;
    end
end

function BFGS = BFGS(B,p,q)
    theta = 1;
    if p'*q < 0.2*p'*B*p
        theta = (0.8*p'*B*p)/(p'*B*p - p'*q);
    end
    r = theta*q + (1-theta)*B*p;
    BFGS = B + (r*r')/(p'*r) - (B*p)*(B*p)'/(p'*B*p);
end