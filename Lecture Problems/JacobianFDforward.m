function J = JacobianFDforward(fun,x)
   n = length(x);
   f = fun(x);
   m = length(f);
   J = zeros(m,n);
   for i=1:n
       h = sqrt(eps);
       x1 = x;
       x1(i) = x(i)+h;
       f1 = feval(fun,x1);
       J(:,i) = (f1-f)/h;
   end
end