function [z,error,A,A1] = smoothing_solver(x,x0,L,SNR,ss)

    [N,n2] = size(x);

    if exist('Params')                == 0,  Params.n1          = N;                end
    if isfield(Params, 'n2')          == 0,  Params.n2          = n2;               end

    if isfield(Params, 'L')           == 0,  Params.L           = L;                end
    if isfield(Params, 'T')           == 0,  Params.T           = 800;              end

    if isfield(Params, 'y1')          == 0,  Params.y1          = 0.1;              end
    if isfield(Params, 'u0')          == 0,  Params.u0          = 65;               end
    if isfield(Params, 'y')           == 0,  Params.y           = 0.1;              end

    if isfield(Params, 'mu')          == 0,  Params.mu          = 0.6;              end
    if isfield(Params, 'B')           == 0,  Params.B           = N*n2;             end

    m           = N*n2*ceil(N*n2/Params.L);
    Params.m    = m;

    display(Params)

    % Make linear operators;
    A  = @(I) fftshift(fft(FROG_signal(I,L,N)),1);       % propagation operator
    A1 = @(J) fftshift(fft(FROG_signal(J,1,N)),1);       % propagation operator

    %% Make signal and data (noiseless)
    y      = abs(A(x)).^2;
        
    if SNR
        Ynoisy = awgn(y,SNR,'measured',ss);
    else
        Ynoisy = y;
    end
    
    fprintf('simulated snr = %f\n',snr(y,y-Ynoisy));
    
    tic
    [z,error] = solver(Ynoisy,x,x0,Params, A,1e-6);
    toc
    z      = best_sol(z, x);
    z  = exp(-1i * angle(trace(x' * z)))*z;
end
