function y = aux(x,L,N,fi)
    time = (-N/2:N/2-1).';
    time = time*1e-15;
    dt   = time(2)-time(1);
    F    = (-N/2:N/2-1).';
    F    =  F/dt/N;
    ind  = 1:N; 
    D    = time(ind).';

    y    = zeros(N,length(fi));
    X    = fft(x);

    for ik=1:size(y,2)
        y(:,ik) = ifft( X.*exp(1i*2*pi*D((fi(ik)-1)*L+1)*F) ) + ifft( X.*exp(-1i*2*pi*D((fi(ik)-1)*L+1)*F).*exp(1i*2*pi*D((fi(ik)-1)*L+1)*F(fi(ik))) );
    end
end