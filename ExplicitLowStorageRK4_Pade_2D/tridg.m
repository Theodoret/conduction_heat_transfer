function x = tridg(a,b,c,d,N)
% Thomas algorithm for tridiagonal matrix solution 
% https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
n = N;
   x=d;
    % forward elimination
    for i = 2:n
        w = a(i)/b(i-1);
        b(i) = b(i)-w*c(i-1);
        d(i) = d(i)-w*d(i-1);     
    end
    
    % backward substitution
    x(n) = d(n)/b(n);
    for i = (n-1):(-1):1
        x(i) = (d(i)-c(i)*x(i+1))/b(i);
    end
    
end