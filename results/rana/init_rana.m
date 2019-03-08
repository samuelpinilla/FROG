function x0 = init_rana(Y)
    %% initialization
    N       = size(Y,1);
    g       = sum(Y,2);
    
    s = sqrt(fftshift(ifft(g))*sqrt(N));
    s = s/norm(s,'fro');
    
    %% variables
    alpha = 0.09;
    beta  = 0.425;
    gamma = 1;
    
    for t=2:N
        
        term1p = alpha*abs(s(t)-s(t-1))^2;
        term1n = alpha*abs(-s(t)-s(t-1))^2;
        
        if t-2 <= 0
            term2p = 0;
            term2n = 0;
        else
            term2p = beta*abs(s(t)-s(t-1)-(s(t-1)-s(t-2)))^2;
            term2n = beta*abs(-s(t)-s(t-1)-(s(t-1)-s(t-2)))^2;
        end
        
        if t-3 <= 0
            term3p = 0;
            term3n = 0;
        else
            term3p = gamma*abs(s(t)-s(t-1)-(s(t-1)-s(t-2)) - (s(t-1)-s(t-2)-(s(t-2)-s(t-3))))^2;
            term3n = gamma*abs(-s(t)-s(t-1)-(s(t-1)-s(t-2)) - (s(t-1)-s(t-2)-(s(t-2)-s(t-3))))^2;
        end
        
        epsilon1 = term1p + term2p + term3p;
        epsilon2 = term1n + term2n + term3n;
        
        if epsilon1 > epsilon2
            s(t) = -s(t);
        end
    end
    x0 = fft(fftshift(s))/sqrt(N);
end