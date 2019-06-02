clc;
clear;
close all;

%% variables

load('pulse_set.mat');

x     = pulse_set(2,:).';
[N,~] = size(x);
L     = 2;

A = @(I) fftshift(fft(FROG_signal(I,L,N)),1);       % propagation operator
B = @(I) fftshift(fft(FROG_signal(I,1,N)),1);       % propagation operator

%% Make signal and data (noiseless)

y  = abs(A(x)).^2;
Y  = zeros(N,N);

SNR = 0;

if SNR > 0
    aux  = awgn(y,SNR,'measured');
    for ii=1:size(y,1)
        Y(ii,:) = interp1(1:L:N, aux(ii,:), 1:N, 'pchip');
    end
    Y = sqrt(Y.*(Y>=0));
else
    for ii=1:size(y,1)
        Y(ii,:) = interp1(1:L:N, y(ii,:), 1:N, 'pchip');
    end
    Y = sqrt(Y.*(Y>=0));
end

z  = init_pg(Y);

error   = norm(Y-abs(B(z)),'fro')/norm(Y,'fro');
fprintf('Iter = 0  Error = %f \n', error);

z = z.';
[z1,z,er] = ispecshg(Y,z,conj(z),B,300);

z  = best_sol(z.', x);

%% plots
T = length(er);

I_in        = quickscale(abs(x).^2);
phase_in    = angle(x);        phase_in(I_in<max(I_in)*5e-3) = NaN;
phase_in    = unwrap(phase_in);

I_ret_s     = quickscale(abs(z).^2);
phase_ret   = angle(z);     phase_ret(I_in<max(I_in)*5e-3) = NaN;
phase_ret_s = unwrap(phase_ret);

figure;
subplot(2,2,1);imagesc(Y.^2),title('Simulated trace');
subplot(2,2,2);imagesc(abs(B(z)).^2),title('Reconstructed trace');
subplot(2,2,3);plotyy(time,I_in,time,I_ret_s),title('Reconstructed amplitude');
subplot(2,2,4);plotyy(time*1e15,phase_in,time*1e15,phase_ret_s),title('Reconstructed phase');
