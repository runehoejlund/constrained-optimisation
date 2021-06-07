function [x, lambda, f, k] = SQPLocalSolver(fun, constraints, x0, lambda0, B0)
    B = B0;
    k = 1;
    x(:,k) = x0;
    lambda(:,k) = lambda0;
    [f(k), df(:,k)] = fun(x(:,k));
    [c(:, k), dc] = constraints(x(:,k));
    
    % Parameters for convergence criteria
    epsilon = 10^(-6);
    
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
        lambda(:,k+1) = lambdas.ineqlin;
        x(:,k+1) = x(:,k) + p;
        
        % Compute gradient of Lagrangian in x_k
        dLk = df(:,k) - dc*lambda(:,k+1);
        
        % Evaluate objective and constraints in x_(k+1)
        [f(k+1), df(:,k+1)] = fun(x(:,k+1));
        [c(:, k+1), dc] = constraints(x(:,k+1));
        
        % Compute gradient of Lagrangian in x_(k+1)
        dLkp1 = df(:,k+1) - dc*lambda(:,k+1);
        q = dLkp1 - dLk;
        
        % Damped BFGS Updating
        theta = 1;
        if p'*q < 0.2*p'*B*p
            theta = (0.8*p'*B*p)/(p'*B*p - p'*q);
        end
        r = theta*q + (1-theta)*B*p;
        B = B + (r*r')/(p'*r) - (B*p)*(B*p)'/(p'*B*p);
        
        k = k + 1;
    end
end