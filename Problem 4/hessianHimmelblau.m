function h = hessianHimmelblau(x,lambda)
    % Hessian of f 
    hessf = [8*x(1)^2 + 4*(x(1)^2+x(2)-11) + 2, 4*x(1)+4;
            4+4*x(2),                          2+8*x(2)^2+4*(x(1)+x(2)^2-7)];
    % Hessian of c
    hessc1 = [2,    0;
              0,    0];
    hessc2 = zeros(2);
    h = hessf + lambda.ineqnonlin(1)*hessc1  + lambda.ineqnonlin(2)*hessc2;
end