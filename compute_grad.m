function grad = compute_grad(z, u, y, Params,i)

    L      = Params.L;
    N      = Params.n1;
    batch  = Params.B;
    
    r      = batch/N;
    rc     = ceil(i/N);
    fi     = rc:rc+r-1;
    
    zf     = aux(z,L,N,fi);
    
    ysub = y(:,fi);
    
    Az   = fftshift(fft(FROG_signal_B(z,L,N,fi)),1);
    yz   = sqrt(abs(Az).^2+u^2);
    aa   = max(max(yz),1);
    grad = 4*sum(conj(zf).*ifft(fftshift(((yz-ysub)./yz).*Az,1)),2)*N/(aa*batch);
end
