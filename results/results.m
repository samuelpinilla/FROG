clc;
clear;
close all;
clear classes;

addpath(genpath('smoothing'));
addpath(genpath('ptych'));
addpath(genpath('pcpg'));

load('pulse_set.mat');

%% Fig. 1
pp = [8,1];
SNR = 0;
figure;

for t = 1:length(pp)
    x = pulse_set(pp(t),:).';
    [z_s,~,~,~] = smoothing_solver(x,[],1,SNR);
    [z_p,~]     = pytch_solver(x,[],1,SNR);

    subplot(length(pp),2,length(pp)*(t-1)+1),plot(time*1e15,abs(x),time*1e15,abs(z_s),time*1e15,abs(z_p)),...
        title('Reconstructed magnitud'), xlabel('Time [fsec]','FontSize',16); ylabel('Amplitude','FontSize',16);
    subplot(length(pp),2,length(pp)*(t-1)+2),plot(time*1e15,unwrap(angle(x)),time*1e15,unwrap(angle(z_s)),...
        time*1e15,unwrap(angle(z_p)))), title('Reconstructed phase'),xlabel('Time [fsec]','FontSize',16); ylabel('Phase','FontSize',16);
end

%% Fig. 2
L   = [2,4,8];
SNR = 20;
figure;
for t=1:length(L)
    x = pulse_set(2,:).';
    [z_s,~,~,A1] = smoothing_solver(x,[],L(t),SNR);
    [z_p,~]      = pytch_solver(x,[],L(t),SNR);

    subplot(length(L),3,length(L)*(t-1)+1),imagesc(abs(A1(z_p))),title('Pytch'),xlabel('Time [fsec]','FontSize',16);
    subplot(length(L),3,length(L)*(t-1)+3),imagesc(abs(A1(z_s))),title('Proposed'),xlabel('Time [fsec]','FontSize',16);
end
%% Fig. 3
L   = [2,4,8];
SNR = 20;
figure;
for t=1:length(L)
    x = pulse_set(2,:).';
    [z_s,~,~,A1] = smoothing_solver(x,[],L(t),SNR);
    [z_p,~]      = pytch_solver(x,[],L(t),SNR);

    subplot(length(L),2,2*(t-1)+1),plot(time*1e15,abs(x),time*1e15,abs(z_s),time*1e15,abs(z_p)),...
        title('Reconstructed magnitud'),xlabel('Time [fsec]','FontSize',16); ylabel('Amplitude','FontSize',16);
    subplot(length(L),2,2*(t-1)+2),plot(time*1e15,unwrap(angle(x)),time*1e15,unwrap(angle(z_s)),...
        time*1e15,unwrap(angle(z_p)))), title('Reconstructed phase'),xlabel('Time [fsec]','FontSize',16); ylabel('Phase','FontSize',16);
end

%% Fig. 4
L   = 1:1:6;
SNR = 0;
prob_s = zeros(12,6);
prob_p = zeros(12,6);
prob_pg = zeros(12,6);
c = 0;

for t=1:length(L)
    x = pulse_set(2,:).';
    for del=1.2:-0.1:0.1
        c = c + 1;
        for it=1:100
            x0  = x + del*randsrc(1,1,[-1,1]);
            [~,error_s,~,~] = smoothing_solver(x,x0,L(t),SNR);
            [~,error_p]     = pytch_solver(x,x0,L(t),SNR);
            [~,error_pg]    = PCPG_sol(x,x0,L(t),SNR);

            if min(error_s(error_s>0))<=1e-6
                prob_s(c,t) = prob_s(c,t) + 1;
            end
            if min(error_p(error_p>0))<=1e-6
                prob_p(c,t) = prob_p(c,t) + 1;
            end
            if min(error_p(error_pg>0))<=1e-6
                prob_pg(c,t) = prob_pg(c,t) + 1;
            end
        end
    end
    c = 0;
end

figure;
subplot(2,2,2*(t-1)+1),imagesc(prob_s/100),colormap gray,title('Proposed'),...
    xlabel('L','FontSize',16); ylabel('delta','FontSize',16);
subplot(2,2,2*(t-1)+2),imagesc(prob_p/100),colormap gray,title('Pytch'),...
    xlabel('L','FontSize',16); ylabel('delta','FontSize',16);
subplot(2,2,2*(t-1)+4),imagesc(prob_pg/100),colormap gray,title('PCPG'),...
    xlabel('L','FontSize',16); ylabel('delta','FontSize',16);

%% Fig. 5
L   = 1:1:8;
SNR = 0;
ferror_s = zeros(length(L),1);
ferror_p = zeros(length(L),1);

for t=1:length(L)
    x = pulse_set(2,:).';
    for it=1:100
        [~,~,error_s,error_p] = smoothing_init(x,L(t),SNR);
        ferror_s(t) = ferror_s(t) + error_s;
        ferror_p(t) = ferror_p(t) + error_p;
    end
    fprintf('L = %d\n',t);
end
figure;
plot(L,ferror_p/100,L,ferror_s/100,L,ferror_r/100),xlabel('L','FontSize',16); ylabel('relative error','FontSize',16);

%% Fig. 6
L   = 4;
SNR = 0;
x   = pulse_set(2,:).';
[z_s,z_p,~,~] = smoothing_init(x,L,SNR);

[z_s,~,~,~] = smoothing_solver(x,z_s,L,SNR);
[z_p,~,~,~] = smoothing_solver(x,z_p,L,SNR);
[z_r,~,~,~] = smoothing_solver(x,z_r,L,SNR);

subplot(1,2,1),plot(time*1e15,abs(x),time*1e15,abs(z_s),time*1e15,abs(z_p)),title('Reconstructed magnitud'),...
    xlabel('Time [fsec]','FontSize',16); ylabel('Amplitude','FontSize',16);
subplot(1,2,2),plot(time*1e15,unwrap(angle(x)),time*1e15,unwrap(angle(z_s)), time*1e15,unwrap(angle(z_p)))),...
    title('Reconstructed phase'),xlabel('Time [fsec]','FontSize',16); ylabel('Phase','FontSize',16);

%% Fig. 7

L   = 1:1:8;
SNR = 8:4:20;
ferror_s = zeros(length(SNR),length(L));

for t=1:length(L)
    x = pulse_set(2,:).';
    for ss=1:length(SNR)
        for it=1:100
            [~,~,error_s,~] = smoothing_init(x,L(t),SNR(ss));
            ferror_s(ss,t)  = ferror_s(ss,t) + error_s;
        end
    end
    fprintf('L = %d SNR = %f\n',t,SNR(ss));
end
figure;
plot(ferror_s/100),xlabel('L','FontSize',16); ylabel('relative error','FontSize',16);
