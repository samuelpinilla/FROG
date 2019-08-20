%% Generate time and frequency
N = 256/2;
Np = 100;
T = 100;
time = linspace(-T,T,N);
time = time*1e-15;
dt = time(2)-time(1);
F = (-N/2:N/2-1);
F =  fftshift( F/dt/N );
SNR = 20;
% Load pulse bank, 100 random pulses
load('pulse_set.mat');
alg = 1; % use 1 for pulse reconstruction without any prior information about the pulse
               % use 2 for pulse reconstruction with known power spectrum
LPF_flag = 0; % 1 for applying Low Pass Filter. 0 without LPF.
%% Make SHG FROG trace
    pulse = pulse_set(32, :);% chooses pulse #2. You can choose another pulse from the bank or use your own pulse here
    ind = 1:N;
    Ns = numel(ind); 
    D = time(ind);
    I = zeros(N, Ns);
    for ik=1:Ns
        P = ifft( fft(pulse).*exp(1i*2*pi*D(ik)*F) );
%         figure(2);plot(time, abs(P), time, abs(pulse), time, abs(P.*pulse) );waitforbuttonpress;
        I(:, ik) = fftshift( abs(fft(P.*pulse)/N) );
    end
Inoisy = awgn(I, SNR, 'measured');

if LPF_flag
    Fsupp = (abs(F)<2.5e13)';
else
    Fsupp = logical(ones(size(F)))';
end
LPF = kron(ones(1,N), fftshift(Fsupp));
InoisyLPF = Inoisy.*LPF;
etta = sum(Fsupp)/N;

ref = conj(pulse');
Terror = inf;
%%  Run algorithm
for ni=1:1   % number of random realizations
    if alg==1
        [Obj, error, Irec] = ePIE_fun_FROG(InoisyLPF, D, 200, Fsupp, F', time, 1e-6);
    else
        [Obj, error, Irec] = ePIE_fun_FROG_sp(InoisyLPF, D, 200, Fsupp, F', time, 1e-3, abs( fft(ref)/1 ));
    end  
            flag = error(find(error, 1, 'last'))<1e-3;
            if flag
                ObjB = best_sol(Obj, ref);
                Aerror = acos(abs(ObjB'*ref)/sqrt( (ObjB'*ObjB)*(ref'*ref) ));
                break
            end
            if Terror > error(find(error, 1, 'last'))
                Terror = error(find(error, 1, 'last'));
                ObjB = best_sol(Obj, conj(pulse'));
                Aerror = acos(abs(ObjB'*ref)/sqrt( (ObjB'*ObjB)*(ref'*ref) ));
            end
 end
%% Plot figures
delta = acos(abs(ObjB'*ref)/sqrt( (ObjB'*ObjB)*(ref'*ref) ));
figure(2)

subplot(3, 2, 1); imagesc(time*1e15, fftshift(F)*1e-15, sqrt(abs(I)))
xlabel('Time [fsec]','FontSize',16); ylabel('Freq.[THz]','FontSize',16);
title('Simulated trace')

subplot(3, 2, 2); imagesc(time*1e15, fftshift(F)*1e-15, sqrt(abs(Inoisy)))
xlabel('Time [fsec]','FontSize',16); ylabel('Freq.[THz]','FontSize',16);
title('Simulated noisy trace')

subplot(3,2,3)
plot(time*1e15, abs(ObjB))
xlabel('Time [fsec]','FontSize',16); ylabel('Amplitude [a.u.]','FontSize',16);
hold all
plot(time*1e15, abs(ref))
xlabel('Time [fsec]','FontSize',16); ylabel('Amplitude [a.u.]','FontSize',16);
subplot(3,2,4)
plot(time*1e15, unwrap(   angle(ObjB)   )/pi)
xlabel('Time [fsec]','FontSize',16); ylabel('Phase [pi]','FontSize',16);
hold all
plot(time*1e15, unwrap(   angle(ref)   )/pi)
xlabel('Time [fsec]','FontSize',16); ylabel('Phase [pi]','FontSize',16);

Icut = zeros(N,N);
indF = find(fftshift(Fsupp));
Icut(indF, ind) = I(indF,:);
subplot(3,2,5)
imagesc(time*1e15, fftshift(F)*1e-15, sqrt(abs(Icut)))
xlabel('Time [fsec]','FontSize',16); ylabel('Freq.[THz]','FontSize',16);
title('Used trace')

subplot(3, 2, 6); imagesc(time*1e15, fftshift(F)*1e-15, sqrt(abs(Irec)))
xlabel('Time [fsec]','FontSize',16); ylabel('Freq.[THz]','FontSize',16);
title('Reconstructed trace')
suptitle( sprintf('delta = %.3f', delta) )
