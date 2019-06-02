clc;
clear;
close all;
clear classes;

%% variables
load('pulse_set.mat');

x       = pulse_set(8,:).';
[n1,n2] = size(x);

if exist('Params')                == 0,  Params.n1          = n1;               end
if isfield(Params, 'n2')          == 0,  Params.n2          = n2;               end

if isfield(Params, 'L')           == 0,  Params.L           = 1;                end
if isfield(Params, 'T')           == 0,  Params.T           = 500;              end

if isfield(Params, 'y1')          == 0,  Params.y1          = 0.3;              end
if isfield(Params, 'u0')          == 0,  Params.u0          = 80;               end
if isfield(Params, 'y')           == 0,  Params.y           = 0.5;              end

if isfield(Params, 'mu')          == 0,  Params.mu          = 0.6;              end
if isfield(Params, 'B')           == 0,  Params.B           = (n1*n2);          end

L           = Params.L;
m           = n1*n2*ceil(n1*n2/Params.L);
Params.m    = m;

display(Params)

% Make linear operators;
A = @(I) fftshift(fft(FROG_signal(I,L,N)),1);       % propagation operator
B = @(I) fftshift(fft(FROG_signal(I,1,N)),1);       % propagation operator

%% Make signal and data (noiseless)
y      = abs(A(x)).^2;
SNR    = 0;

if SNR
    Ynoisy = awgn(y,SNR,'measured');
else
    Ynoisy = y;
end

tic
[z,error,normgrad] = solver(Ynoisy,x,Params, A,1e-6);
toc

z      = best_sol(z, x);
delta  = min(error);

I_in        = quickscale(abs(x).^2);       phase_in    = angle(x);
phase_in(I_in<max(I_in)*5e-3) = NaN;         
phase_in    = unwrap(phase_in);   

I_ret_s     = quickscale(abs(z).^2);     
phase_ret   = angle(z);
phase_ret(I_in<max(I_in)*5e-3) = NaN;      
phase_ret_s = unwrap(phase_ret);

%% plots

ff = figure;
axes1 = axes('Parent',ff,...
    'Position',[0.13 0.720906497767073 0.334659090909091 0.135967246868178]);
hold(axes1,'on');

subplot(2,2,1);imagesc(Ynoisy),title('Simulated trace');
subplot(2,2,2);imagesc(abs(B(z)).^2),title('Reconstructed trace');

subplot(2,2,3);[~,H1,H2] = plotyy(time,I_ret_s,time,I_in),title('Reconstructed amplitude');
set(H1,'LineWidth',2,'Color',[1 0.800000011920929 0.400000005960464]);
set(H2,'LineWidth',1,'LineStyle','-','Color',[1 0 0]);

subplot(2,2,4);[~,H1,H2] = plotyy(time*1e15,phase_ret_s,time*1e15,phase_in),title('Reconstructed phase');
set(H1,'LineWidth',2,'Color',[1 0.800000011920929 0.400000005960464]);
set(H2,'LineWidth',1,'LineStyle','-','Color',[1 0 0]);