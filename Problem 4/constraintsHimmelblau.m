function [c, dc] = constraintsHimmelblau(x)
    xl = -5;
    xu = 5;
    g = [(x(1)+2).^2 - x(2);
        -4*x(1)+10*x(2)];
    dg = [2*(x(1)+2)      -4;
            -1            10];
    c = [g; x(1)-xl; x(2)-xl; x(1) + xu; x(2) + xu];
    dc = [dg eye(2) -eye(2)];
end