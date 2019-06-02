function y = FROG_signal(x,L,N)
%% variables
    time = (-N/2:N/2-1).';
    time = time*1e-15;
    dt   = time(2)-time(1);
    F    = (-N/2:N/2-1).';
    F    =  fftshift( F/dt/N );
    ind  = 1:N; 
    D    = time(ind).';

    y    = zeros(N,ceil(N/L));
    X    = fft(x);

    for ik=1:size(y,2)
        P       = ifft( X.*exp(1i*2*pi*D((ik-1)*L+1)*F) );
        y(:,ik) = P.*x;
    end
end