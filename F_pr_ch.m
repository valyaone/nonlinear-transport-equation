function f_vector = F_pr_ch(U,t,h)
    
    Uleft = @(t) exp(-t);

    N = length( U );
    
    f_vector(1) = -1/(h)*( U(1)^2 - U(1)*Uleft(t) ) + exp( U(1)^2 );

    for n=2:N
        f_vector(n) = -1/(h)*(U(n)^2 - U(n)*U(n-1)) + exp(U(n)^2);
    end
    
end