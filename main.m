clc;
clear;
close all;

%% variables

load('pulse_set.mat');

x       = pulse_set(2,:).';
[n1,n2] = size(x);

if exist('Params')                == 0,  Params.n1          = n1;               end
if isfield(Params, 'n2')          == 0,  Params.n2          = n2;               end

if isfield(Params, 'L')           == 0,  Params.L           = 2;                end
if isfield(Params, 'T')           == 0,  Params.T           = 300;              end

if isfield(Params, 'y1')          == 0,  Params.y1          = 0.3;              end
if isfield(Params, 'u0')          == 0,  Params.u0          = 80;               end
if isfield(Params, 'y')           == 0,  Params.y           = 0.5;              end

if isfield(Params, 'mu')          == 0,  Params.mu          = 0.8;              end
if isfield(Params, 'B')           == 0,  Params.B           = n1*n2;            end

L           = Params.L;
m           = n1*n2*ceil(n1*n2/Params.L);
Params.m    = m;

display(Params)

% Make linear operators;
A = @(I) fftshift(fft(FROG_signal(I,L,N)),1);       % propagation operator

%% Make signal and data (noiseless)
y      = abs(A(x));

SNR    = 0;

if SNR
    Ynoisy = awgn(y,SNR,'measured');
else
    Ynoisy = y;
end

tic
[z,error] = solver(Ynoisy,x,Params, A,1e-6);
toc


z      = best_sol(z, x);
delta  = min(error);


%% plots
T = length(error)-1;

figure;
subplot(4,2,1);imagesc(y),title('Simulated trace');
subplot(4,2,2);imagesc(Ynoisy),title('Simulated noisy trace');
subplot(4,2,3);plot(time,abs(x),time,abs(z)),title('Reconstructed amplitude');
subplot(4,2,4);plot(time*1e15,unwrap(angle(x))/pi,time*1e15,unwrap(angle(z))/pi),title('Reconstructed phase');
subplot(4,2,5);imagesc(Ynoisy),title('Used trace');
subplot(4,2,6);imagesc(abs(A(z))),title('Reconstructed trace'),suptitle( sprintf('delta = %.3f', delta) );
subplot(4,2,7);semilogy(0:T,error),xlabel('Iteration'), ylabel('Relative error (log10)'), ...
title('Relative error vs. iteration count')