function [f, df, d2f] = Differentiate(fun,x)
    f = fun(x);
    n = length(x); m = length(f);
    df = zeros(n,m); d2f = zeros(n,n,m); 
    E = eye(n);
    for i=1:n
       h = sqrt(eps);
       dxi = x + h*E(:,i);
       df(i,:) = (fun(dxi) - f)/h;
       for j=1:n
           h = sqrt(h);
           dxi = x + h*E(:,i);
           if i==j
               dxj = x - h*E(:,j);
               d2f(i,j,:) = (fun(dxi)+fun(dxj)-2*f)/(h^2);
               continue;
           end
           dxj = x + h*E(:,j);
           dxij = dxi + dxj - x;
           d2f(i,j,:) = (fun(dxij) - fun(dxi) - fun(dxj) + f)/(h^2);
       end
    end
end