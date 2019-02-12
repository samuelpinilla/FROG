function g = init_pg(Y, L,A)
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
    end
end