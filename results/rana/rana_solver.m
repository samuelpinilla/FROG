function [z,er] = rana_solver(x,x0,L,SNR,filter)
%% variables
[N,~] = size(x);

A  = @(I) fftshift(fft(FROG_signal(I,1,N)),1);       % propagation operator

%% Make signal and data (noiseless)

y     = abs(A(x));
ind_L = 1:L:N;

if SNR >0
    aux = awgn(y,SNR,'measured');
    Ynoisy = zeros(size(y));
    Ynoisy(:,ind_L) = aux(:,ind_L);
else
    Ynoisy = y;
end

if isempty(x0)
    z  = init_rana(Ynoisy.^2,filter);
else
    z = x0;
end

error   = norm(y-abs(A(z)),'fro')/norm(y,'fro');
fprintf('Iter = 0  Error = %f \n', error);

z = z.';

tic
[~,z,er] = ispecshg(Ynoisy,z,z,A,300);
toc

z  = best_sol(z.', x);

%% plots
% T = length(er);
% 
% figure;
% subplot(2,2,1);imagesc(y),title('Simulated trace');
% subplot(2,2,2);imagesc(abs(A(z))),title('Reconstructed trace');
% subplot(2,2,3);plot(time,abs(x),time,abs(z)),title('Reconstructed amplitude');
% subplot(2,2,4);plot(time*1e15,unwrap(angle(x))/pi,time*1e15,unwrap(angle(z))/pi),title('Reconstructed phase');

end
