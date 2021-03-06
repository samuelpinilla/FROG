function x0 = init(Y,L)
%% initialization
    N       = size(Y,1);
    Z       = fft(fftshift(Y,1))./N;
    Yhat    = zeros(N,N);
    
    if L > 1
        for ii=1:size(Y,1)
            Yhat(ii,:) = interp1(1:L:N, Z(ii,:), 1:N, 'pchip');
        end
        Y1      = sqrt(abs(real(fftshift(ifft(Yhat),1))));
        
        g       = sum(Y1,1)'/sqrt(sum(sum(abs( sum(Y1,1) ).^2)));
        g       = g.*exp(-1i*randsrc(N,1,[0,1],275));
    else
        g       = sum(Y,1)'/sqrt(sum(sum(abs( sum(Y,1) ).^2)));
        
        pp      = exp(-1i*randsrc(N,1,[0,1],275));
        g       = g.*pp;
        Yhat    = Z;
    end
    I       = eye(N);
        
    G       = zeros(size(Yhat));
    
    time    = (-N/2:N/2-1).';
    time    = time*1e-15;
    dt      = time(2)-time(1);
    F       = (-N/2:N/2-1).';
    F       = F/dt/N;
    ind     = 1:N; 
    D       = time(ind).';

    X       = fft(g)/sqrt(N);
    X0      = zeros(N,N);
    
    for iter=1:2
        for il=1:N
            gl     = g.*conj(ifft( X.*exp(1i*2*pi*D(il)*F ))*sqrt(N));
            for ik=1:N
                P1      = ifft( X.*exp(1i*2*pi*D(ik)*F) )*sqrt(N);
                P2      = ifft( X.*exp(1i*2*pi*D(mod(ik+il-2,N)+1)*F ))*sqrt(N);
                G(ik,:) = P1.*conj(P2);
            end
            
            zl       = Yhat(il,:).';
            X0(:,il) = (G'*G + 0.5*I)\( 0.5*gl + G'*zl);
        end
        xd          = convert2diag(X0);
        [eigvec, ~] = eigs(xd,1);
        g           = eigvec;
        X           = fft(g)/sqrt(N);
        
    end
    d    = diag(xd);
    x0   = power(sum(d.*(d>0)),1/4)*g;
end