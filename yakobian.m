function f_vector = yakobian(U,t,h)
    
    Uleft = @(t) exp(-t);

    N = length( U );
    
    f_vector(1,1) = -1/h*( 2*U(1) - Uleft(t) ) + 2*U(1) *exp( U(1)^2 );
        
    for n=2:N
        f_vector(n,n-1) = U(n)/h;
        f_vector(n,n) = -1/h*(2*U(n) - U(n-1)) + 2*U(n) *exp(U(n-1)^2);
    end
    
end