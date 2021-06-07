function [fminc, ceq, fmindc, dceq] = fminconstraintsHimmelblau(x,p)
    [c, dc] = constraintsHimmelblau(x);
    fminc = -c;
    fmindc = -dc;
    ceq = [];
    dceq = [];
end