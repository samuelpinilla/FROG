function [z_s,z_p,error_s,error_p] = smoothing_init(x,L,SNR)

    [N,n2] = size(x);

    if exist('Params')                == 0,  Params.n1          = N;                end
    if isfield(Params, 'n2')          == 0,  Params.n2          = n2;               end
    if isfield(Params, 'L')           == 0,  Params.L           = L;                end
    if isfield(Params, 'B')           == 0,  Params.B           = N*n2;             end

    m           = N*n2*ceil(N*n2/Params.L);
    Params.m    = m;

    % Make linear operators;
    A  = @(I) fftshift(fft(FROG_signal(I,L,N)),1);       % propagation operator

    %% Make signal and data (noiseless)
    y      = abs(A(x));
    
    if SNR>0
       Ynoisy = awgn(y,SNR,'measured');
    else
       Ynoisy = y;
    end

    [z_s,z_p]  = init(Ynoisy,Params.L,A);
    error_s    = norm(Ynoisy-abs(A(z_s)),'fro')/norm(Ynoisy,'fro');
    z_s        = best_sol(z_s, x);
    
    error_p    = norm(Ynoisy-abs(A(z_p)),'fro')/norm(Ynoisy,'fro');
    z_p        = best_sol(z_p, x);

end