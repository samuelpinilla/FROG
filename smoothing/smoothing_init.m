function [z_s,error_s] = smoothing_init(x,L,SNR,ss)

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
    y      = abs(A(x)).^2;
    
    rng(ss.Seed,ss.Type); % for reproducibility
    switch SNR
        case 0
            Ynoisy = y;
        case 20
            alpha  = 0.1;
            m      = randn(size(y));
            hh     = 1+m*alpha;
            Ynoisy = y.*hh;
            fprintf('simulated snr = %f\n',snr(y,y.*m.*alpha));
        case 16
            alpha  = 0.15;
            m      = randn(size(y));
            hh     = 1+m*alpha;
            Ynoisy = y.*hh;
            fprintf('simulated snr = %f\n',snr(y,y.*m.*alpha));
        case 12
            alpha  = 0.25;
            m      = randn(size(y));
            hh     = 1+m*alpha;
            Ynoisy = y.*hh;
            fprintf('simulated snr = %f\n',snr(y,y.*m.*alpha));
        case 8
            alpha  = 0.35;
            m      = randn(size(y));
            hh     = 1+m*alpha;
            Ynoisy = y.*hh;
            fprintf('simulated snr = %f\n',snr(y,y.*m.*alpha));
    end


    z_s        = init(Ynoisy,Params.L);
    error_s    = norm(sqrt(Ynoisy)-abs(A(z_s)),'fro')/norm(sqrt(Ynoisy),'fro');
    z_s        = best_sol(z_s, x);

end