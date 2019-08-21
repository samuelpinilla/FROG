clc;
clear;
close all;

addpath(genpath('smoothing'));
addpath(genpath('ptych'));

load('pulse_set.mat');
load('seed.mat'); % for reproducibility

%% Fig. 1
pp = [8,1];
SNR = 0;
figure;

for t = 1:length(pp)
    x = pulse_set(pp(t),:).';
    [z_s,~,~,~]  = smoothing_solver(x,[],1,SNR,ss);
    [z_p,~]      = pytch_solver(x,[],1,SNR,ss);
       
    I_in        = quickscale(abs(x).^2);       
    phase_in    = angle(x);        phase_in(I_in<max(I_in)*5e-3) = NaN;         
    phase_in    = unwrap(phase_in);   
    
    I_ret_s     = quickscale(abs(z_s).^2);     
    phase_ret   = angle(z_s);     phase_ret(I_in<max(I_in)*5e-3) = NaN;      
    phase_ret_s = unwrap(phase_ret);
    
    I_ret_p     = quickscale(abs(z_p).^2);     
    phase_ret   = angle(z_p);     phase_ret(I_in<max(I_in)*5e-3) = NaN;      
    phase_ret_p = unwrap(phase_ret);
    
     
    subplot(length(pp),2,length(pp)*(t-1)+1);
    plot(time*1e15, I_in, time*1e15, I_ret_s, time*1e15, I_ret_p);...
        title('Reconstructed magnitude'); xlabel('Time [fsec]'); ylabel('Amplitude');

    subplot(length(pp),2,length(pp)*(t-1)+2);
    plot(time*1e15, phase_in,time*1e15, phase_ret_s,time*1e15, phase_ret_p);...
        title('Reconstructed phase'); xlabel('Time [fsec]'); ylabel('Amplitude');
end

%% Fig. 2
pp = [8,1];
SNR = 20;
figure;

for t = 1:length(pp)
    x = pulse_set(pp(t),:).';
    [z_s,~,~,~] = smoothing_solver(x,[],1,SNR,ss);
    [z_p,~]     = pytch_solver(x,[],1,SNR,ss);

    I_in        = quickscale(abs(x).^2);       
    phase_in    = angle(x);        phase_in(I_in<max(I_in)*5e-3) = NaN;         
    phase_in    = unwrap(phase_in);  
    
    I_ret_s     = quickscale(abs(z_s).^2);     
    phase_ret   = angle(z_s);     phase_ret(I_in<max(I_in)*5e-3) = NaN;      
    phase_ret_s = unwrap(phase_ret);
    
    I_ret_p     = quickscale(abs(z_p).^2);     
    phase_ret   = angle(z_p);     phase_ret(I_in<max(I_in)*5e-3) = NaN;      
    phase_ret_p = unwrap(phase_ret);
    
    
    subplot(length(pp),2,length(pp)*(t-1)+1);
    plot(time*1e15,I_in,time*1e15,I_ret_s,time*1e15,I_ret_p),...
        title('Reconstructed magnitude'), xlabel('Time [fsec]','FontSize',16); ylabel('Amplitude','FontSize',16);
    subplot(length(pp),2,length(pp)*(t-1)+2);
    plot(time*1e15,phase_in,time*1e15,phase_ret_s,...
        time*1e15,phase_ret_p), title('Reconstructed phase'),xlabel('Time [fsec]','FontSize',16); ylabel('Phase','FontSize',16);
end

%% Fig. 3
L   = [2,4,8];
SNR = 20;
figure;
for t=1:length(L)
    x = pulse_set(2,:).';
    [z_s,~,~,A1] = smoothing_solver(x,[],L(t),SNR,ss);
    [z_p,~]      = pytch_solver(x,[],L(t),SNR,ss);

    I_in        = quickscale(abs(x).^2);       
    phase_in    = angle(x);        phase_in(I_in<max(I_in)*5e-3) = NaN;         
    phase_in    = unwrap(phase_in);   
    
    I_ret_s     = quickscale(abs(z_s).^2);     
    phase_ret   = angle(z_s);     phase_ret(I_in<max(I_in)*5e-3) = NaN;      
    phase_ret_s = unwrap(phase_ret);
    
    I_ret_p     = quickscale(abs(z_p).^2);     
    phase_ret   = angle(z_p);     phase_ret(I_in<max(I_in)*5e-3) = NaN;      
    phase_ret_p = unwrap(phase_ret);
     
    subplot(length(L),2,2*(t-1)+1);
    plot(time*1e15, I_in, time*1e15, I_ret_s, time*1e15, I_ret_p);...
        title('Reconstructed magnitude'); xlabel('Time [fsec]'); ylabel('Amplitude');

    subplot(length(L),2,2*(t-1)+2);
    plot(time*1e15, phase_in,time*1e15, phase_ret_s,time*1e15, phase_ret_p);...
        title('Reconstructed phase'); xlabel('Time [fsec]'); ylabel('Amplitude');
end


%% Fig. 4
L   = [2,4,8];
SNR = 20;
figure;
for t=1:length(L)
    x = pulse_set(2,:).';
    [z_s,~,~,A1]   = smoothing_solver(x,[],L(t),SNR,ss);
    [z_p,~]        = pytch_solver(x,[],L(t),SNR,ss);

    subplot(length(L),2,2*(t-1)+1),imagesc(abs(A1(z_p)).^2),title('Pytch'),xlabel('Time [fsec]','FontSize',16); ylabel('Frequency','FontSize',16);
    subplot(length(L),2,2*(t-1)+2),imagesc(abs(A1(z_s)).^2),title('Proposed'),xlabel('Time [fsec]','FontSize',16); ylabel('Frequency','FontSize',16);
end


%% Fig. 9
L   = [1,2,4,6];
SNR = 0;

iter     = zeros(length(L),1);

time_s   = zeros(length(L),1);

iter_r   = zeros(length(L),1);
prob_sr  = zeros(length(L),1);
prob_s   = zeros(length(L),1);

cr       = zeros(length(L),1);

for ll=1:length(L)
    for t=1:100
        x = pulse_set(t,:).';
        x0 = randn(size(x)) +  1i*randn(size(x));
        
        tic
        [~,error_s,~,~]  = smoothing_solver(x,[],L(ll),SNR,ss);
        
        time_s(ll) = time_s(ll) +  toc;
              
        [~,error_sr,~,~] = smoothing_solver(x,x0,L(ll),SNR,ss);
        
        
        if min(error_s(error_s>0))<=1e-6
            prob_s(ll) = prob_s(ll) + 1;
            iter(ll)   = iter(ll) + length(error_s);
        end
                
        if min(error_sr(error_sr>0))<=1e-6
            cr(ll) = cr(ll) + 1;
            prob_sr(ll) = prob_sr(ll) + 1;
            iter_r(ll)   = iter_r(ll) + length(error_sr);
        end
    end
end

fprintf('number iterations designed = %f\n',iter(1)/100);
fprintf('number iterations random = %f\n',iter(1)/cr(1));
fprintf('time smoothing = %f\n',time_s(1)/100);

figure;
plot(L,prob_s/100),title('Proposed'),...
    xlabel('L','FontSize',16); ylabel('success rate','FontSize',16);


%% Fig. 6
L   = 1:1:8;
SNR = 0;
ferror_s = zeros(length(L),1);
ferror_p = zeros(length(L),1);

for t=1:length(L)
    for it=1:100
        x = pulse_set(it,:).';
        [~,error_s] = smoothing_init(x,L(t),SNR,ss);
        [~,error_p] = pytch_init(x,L(t),SNR,ss);
        
        ferror_s(t) = ferror_s(t) + error_s;
        ferror_p(t) = ferror_p(t) + error_p;
    end
    fprintf('L = %d\n',t);
end
figure;
plot(L,ferror_p/100,L,ferror_s/100),xlabel('L','FontSize',16); ylabel('relative error','FontSize',16);

%% Fig 7.
L   = 4;
SNR = 0;
x   = pulse_set(2,:).';
figure;

[z_s,~] = smoothing_init(x,L,SNR,ss);
[z_p,~] = pytch_init(x,L,SNR,ss);

[z_s,~,~,~]  = smoothing_solver(x,z_s,L,SNR,ss);
[z_p,~,~,~]  = smoothing_solver(x,z_p,L,SNR,ss);

I_in        = quickscale(abs(x).^2);
phase_in    = angle(x);        phase_in(I_in<max(I_in)*5e-3) = NaN;
phase_in    = unwrap(phase_in);

I_ret_s     = quickscale(abs(z_s).^2);
phase_ret   = angle(z_s);     phase_ret(I_in<max(I_in)*5e-3) = NaN;
phase_ret_s = unwrap(phase_ret);

I_ret_p     = quickscale(abs(z_p).^2);
phase_ret   = angle(z_p);     phase_ret(I_in<max(I_in)*5e-3) = NaN;
phase_ret_p = unwrap(phase_ret);
    

subplot(1,2,1);
plot(time*1e15,I_in,time*1e15,I_ret_s,time*1e15,I_ret_p),...
    title('Reconstructed magnitude'), xlabel('Time [fsec]','FontSize',16); ylabel('Amplitude','FontSize',16);
subplot(1,2,2);
plot(time*1e15,phase_in,time*1e15,phase_ret_s,...
    time*1e15,phase_ret_p), title('Reconstructed phase'),xlabel('Time [fsec]','FontSize',16); ylabel('Phase','FontSize',16);



%% Fig. 5
L   = 1:1:6;
SNR = 0;
prob_s = zeros(12,6);
c = 0;

for t=1:length(L)
    x = pulse_set(2,:).';
    for del=1.2:-0.1:0.1
        c = c + 1;
        for it=1:100
            x0  = x + del*randsrc(1,1,[-1,1],275);
            [~,error_s,~,~] = smoothing_solver(x,x0,L(t),SNR,ss);
            

            if min(error_s(error_s>0))<=1e-6
                prob_s(c,t) = prob_s(c,t) + 1;
            end
        end
    end
    c = 0;
end

figure;
imagesc(prob_s/100),colormap gray,title('Proposed'),...
    xlabel('L','FontSize',16); ylabel('delta','FontSize',16);
%% Fig 8.
L   = 1:1:8;
SNR = [0,20,16,12,8];
ferror_s = zeros(length(SNR),length(L));
figure; hold on;

for s=1:length(SNR)
    for t=1:length(L)
        for it=1:100
            x = pulse_set(it,:).';
            [~,error_s] = smoothing_init(x,L(t),SNR(s),ss);
            
            ferror_s(s,t) = ferror_s(s,t) + error_s;
        end
    end
    fprintf('SNR = %d\n',SNR(s));
    plot(L,ferror_s(s,:)/100),xlabel('L','FontSize',16); ylabel('relative error','FontSize',16);
end

