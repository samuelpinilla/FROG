function grad = compute_grad(z, u, y, Params, A,i)
    m      = Params.m;
    L      = Params.L;
    N      = Params.n1;
    batch  = Params.B;
    
    zf       = aux(z,L,N);
    
    ysub = zeros(size(y));
    ysub(i:i+batch-1) = 1;
    
    Az   = A(z);
    yz   = sqrt(abs(Az).^2+u^2);
    grad = 4*sum(conj(zf).*ifft(fftshift(((yz-y)./yz).*Az.*ysub,1)),2)*N/batch;
end
      
