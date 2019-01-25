function x0 = init2(Y, L,A,x)
%% initialization
    N       = size(Y,1);
    Z       = fft(fftshift(Y,1).^2)./N;
    Yhat    = zeros(N,N);
    
    if L > 1
        for ii=1:size(Y,1)
            Yhat(ii,:) = interp1(1:L:N, Z(ii,:), 1:N, 'pchip');
        end
        Y1      = sqrt(abs(real(fftshift(ifft(Yhat),1))));
        
        g       = sum(Y1,1)'/sqrt(sum(sum(abs( sum(Y1,1) ).^2)));
        g       = g.*exp(-1i*randsrc(N,1,[0,1]));
    else
        g       = sum(Y,1)'/sqrt(sum(sum(abs( sum(Y,1) ).^2)));
        g       = g.*exp(-1i*randsrc(N,1,[0,1]));
        Yhat    = Z;
    end
    x0 = g;
    I       = eye(N);
    
    error   = norm(Y-abs(A(g)),'fro')/norm(Y,'fro');
    
    G       = zeros(size(Yhat));
    
    T       = 100;
    time    = linspace(-T,T,N).';
    time    = time*1e-15;
    dt      = time(2)-time(1);
    F       = (-N/2:N/2-1).';
    F       = fftshift( F/dt/N );
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
            X0(:,il) = (G'*G + 0.01*I)\( 0.01*gl + G'*zl);
        end
        xd          = convert2diag(X0);
        [eigvec, ~] = eigs(xd,1);
        g           = eigvec;
        X           = fft(g)/sqrt(N);
        
        error(iter + 1) = norm(Y-abs(A(g)),'fro')/norm(Y,'fro');
    end
    d    = diag(xd);
    x0   = power(sum(d.*(d>0)),1/4)*g;
end