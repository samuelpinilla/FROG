clc;
clear;
close all;

%% variables

load('pulse_set.mat');

x     = pulse_set(8,:).';
[N,~] = size(x);
L     = 1;
ind_L = 1:L:N;

A = @(I) fftshift(fft(FROG_signal(I,L,N)),1);       % propagation operator

%% Make signal and data (noiseless)

y  = abs(A(x));
z  = init_pg(y,L,A);

error   = norm(y-abs(A(z)),'fro')/norm(y,'fro');
fprintf('Iter = 0  Error = %f \n', error);

z = z.';
[z1,z,er] = ispecshg(y,z,z,A,300);

z  = best_sol(z.', x);

%% plots
T = length(er);

figure;
subplot(2,2,1);imagesc(y),title('Simulated trace');
subplot(2,2,2);imagesc(abs(A(z))),title('Reconstructed trace');
subplot(2,2,3);plot(time,abs(x),time,abs(z)),title('Reconstructed amplitude');
subplot(2,2,4);plot(time*1e15,unwrap(angle(x))/pi,time*1e15,unwrap(angle(z))/pi),title('Reconstructed phase');
