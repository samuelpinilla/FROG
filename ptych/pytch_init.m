function [z_p,error_p] = pytch_init(x,L,SNR,ss)

    [N,n2] = size(x);

    if exist('Params')                == 0,  Params.n1          = N;                end
    if isfield(Params, 'n2')          == 0,  Params.n2          = n2;               end
    if isfield(Params, 'L')           == 0,  Params.L           = L;                end
    if isfield(Params, 'B')           == 0,  Params.B           = N*n2;             end

    m           = N*n2*ceil(N*n2/Params.L);
    Params.m    = m;

    % Make linear operators;
    A  = @(I) fftshift(fft(FROG_signal(I,L,N)),1);       % propagation operator
    B  = @(I) fftshift(fft(FROG_signal(I,1,N)),1);       % propagation operator
    %% Make signal and data (noiseless)
    
    y  = abs(A(x)).^2;
    Y  = zeros(N,N);
    
    if SNR
        Ynoisy = awgn(y,SNR,'measured',ss);
        
        for ii=1:size(Ynoisy,1)
            Y(ii,:) = interp1(1:L:N, Ynoisy(ii,:), 1:N, 'pchip');
        end
        Ynoisy = sqrt(Y.*(Y>=0));
    else
        Ynoisy = y;
        for ii=1:size(Ynoisy,1)
            Y(ii,:) = interp1(1:L:N, Ynoisy(ii,:), 1:N, 'pchip');
        end
        Ynoisy = sqrt(Y.*(Y>=0));
    end

    z_p        = sum(Ynoisy,1)'/sqrt(sum(sum(abs( sum(Ynoisy,1) ).^2)));
    error_p    = norm(Ynoisy-abs(B(z_p)),'fro')/norm(Ynoisy,'fro');
    z_p        = best_sol(z_p, x);

end
