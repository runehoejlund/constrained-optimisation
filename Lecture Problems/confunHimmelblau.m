function [c,ceq] = confunHimmelblau(x,p)
    ceq = zeros(0,1);
    c = [-((x(1)+2)^2 - x(2));
              -(-4*x(1) + 10*x(2))]
end