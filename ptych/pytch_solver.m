function [ObjB,error] = pytch_solver(x,x0,L,SNR,ss)
    N = 256/2;
    Np = 100;
    T = 100;
    time = linspace(-T,T,N);
    time = time*1e-15;
    dt = time(2)-time(1);
    F = (-N/2:N/2-1);
    F =  fftshift( F/dt/N );
    % Load pulse bank, 100 random pulses
    alg = 1; % use 1 for pulse reconstruction without any prior information about the pulse
    % use 2 for pulse reconstruction with known power spectrum
    LPF_flag = 0; % 1 for applying Low Pass Filter. 0 without LPF.
    
    %% Make SHG FROG trace
    pulse = x.';% chooses pulse #2. You can choose another pulse from the bank or use your own pulse here
    ind = 1:N;
    Ns = numel(ind);
    D = time(ind);
    I = zeros(N, Ns);
    
    for ik=1:Ns
        P = ifft( fft(pulse).*exp(1i*2*pi*D(ik)*F) );
        %         figure(2);plot(time, abs(P), time, abs(pulse), time, abs(P.*pulse) );waitforbuttonpress;
        I(:, ik) = fftshift( abs(fft(P.*pulse)/N).^2 );
    end
    
    ind_L = 1:L:N;
    
    if SNR > 0       
        aux    = I(:,ind_L);
        Inoisy = awgn(aux,SNR,'measured',ss);
        Inoisy = sqrt(Inoisy.*(Inoisy>=0));
        
    else
        Inoisy = I(:,ind_L);
        Inoisy = sqrt(Inoisy.*(Inoisy>=0));
    end

    if LPF_flag
        Fsupp = (abs(F)<2.5e13)';
    else
        Fsupp = logical(ones(size(F)))';
    end
    
    InoisyLPF = Inoisy;

    ref = conj(pulse');
    Terror = inf;
    %%  Run algorithm
    for ni=1:1   % number of random realizations
        if alg==1
            [Obj, error, Irec] = ePIE_fun_FROG(InoisyLPF,L,x0, D, 300, Fsupp, F', time, 1e-6);
        else
            [Obj, error, Irec] = ePIE_fun_FROG_sp(InoisyLPF, D, 300, Fsupp, F', time, 1e-3, abs( fft(ref)/1 ));
        end
        flag = error(find(error, 1, 'last'))<1e-3;
        if flag
            ObjB = best_sol(Obj, ref);
            ObjB   = exp(-1i * angle(trace(x' * ObjB)))*ObjB;
            Aerror = acos(abs(ObjB'*ref)/sqrt( (ObjB'*ObjB)*(ref'*ref) ));
            break
        end
        if Terror > error(find(error, 1, 'last'))
            Terror = error(find(error, 1, 'last'));
            ObjB = best_sol(Obj, conj(pulse'));
            ObjB   = exp(-1i * angle(trace(x' * ObjB)))*ObjB;
            Aerror = acos(abs(ObjB'*ref)/sqrt( (ObjB'*ObjB)*(ref'*ref) ));
        end
    end
end
