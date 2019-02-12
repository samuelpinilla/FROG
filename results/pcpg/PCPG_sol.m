function [z,er] = PCPG_sol(x,x0,L,SNR)
    %% variables

    [N,~] = size(x);
    ind_L = 1:L:N;

    A = @(I) fftshift(fft(FROG_signal(I,1,N)),1);       % propagation operator

    %% Make signal and data (noiseless)

    y  = abs(A(x));
    
    if SNR > 0
        aux    = awgn(y,SNR,'measured');
        Ynoise = zeros(size(y));
        Ynoise(:,ind_L) = aux(:,ind_L);
    else
        Ynoise = y;
    end
    if isempty(x0)    
        z  = init(Ynoise,L,A);
    else
        z  = x0;
    end

    error   = norm(Ynoise-abs(A(z)),'fro')/norm(Ynoise,'fro');
    fprintf('Iter = 0  Error = %f \n', error);

    z = z.';
    [~,z,er] = ispecshg(Ynoise,z,z,A,300);

    z  = best_sol(z.', x);

    %% plots
%     T = length(er);
%     
%     figure;
%     subplot(2,2,1);imagesc(y),title('Simulated trace');
%     subplot(2,2,2);imagesc(abs(A(z))),title('Reconstructed trace');
%     subplot(2,2,3);plot(time,abs(x),time,abs(z)),title('Reconstructed amplitude');
%     subplot(2,2,4);plot(time*1e15,unwrap(angle(x))/pi,time*1e15,unwrap(angle(z))/pi),title('Reconstructed phase');
end
