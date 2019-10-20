%% ePIE function for SHG FROG 
function [Obj, error, Ir] = ePIE_fun_FROG(I,L,x0, D, iterMax, Fsupp, F, time, STOPc,ind_L,ss)
% reconstructs a pulse function (in time) from a SHG FROG
%   trace by use of the Ptychographic algorithm.
%
%Usage:
%
%   [Obj, error, Ir] = ePIE_fun_FROG_sp(I, D, iterMax, Fsupp, F, time, STOPc, spec)
%
%       Obj      =   Reconstructed pulse field (in time).
%       Ir      =   reconstructed FROG trace.
%       error     =   vector of errors for each iteration
%
%       I       =   Experimental / Simulated SHG FROG Trace
%		D	=	vector of delays that coresponds to trace.
%       iterMax =   Maximum number of iterations allowed (default = 1000).
%       Fsupp = vector of informative frequencies (logical).
%       F = vector of frequencies.
%       time = coresponding vector in time domain.
%       STOPc  =   (Optioanl) Tolerence on the error (default = 1e-5).

%   Set maximum number of iterations
if (~exist('iterMax', 'var')||isempty(iterMax))
    iterMax = 1000;
end

%   Set convergence limit
if (~exist('STOPc', 'var')||isempty(STOPc))
    STOPc = 1e-5;
end


if L>1
    aux = zeros(size(I,1),size(I,1));
    for ii=1:size(I,1)
        aux(ii,:) = interp1(1:L:size(I,1), I(ii,:), 1:size(I,1), 'pchip');
    end
else
    aux = I;
end

[N, K] = size(I);

if isempty(x0)
    Obj = sum(aux,1)'/sqrt(sum(sum(abs( sum(aux,1) ).^2)));
else
    Obj = x0;
end

D = D(ind_L);

del = 1e-3;
del2 = 2e-6;
error = zeros(iterMax, 1);
iter = 1;
Ir = zeros(size(I));

while iter <= iterMax
    rng(5,ss.Type);
    s = randperm(K);
    alpha = abs( 0.2+randn(1,1)/20 );
    for iterK =1:K
        
        temp = sig_shift(Obj, D(s(iterK)), F);
        psi = Obj.*temp;
        psi_n = fft(psi)/N;
        phase = exp(1i*angle(psi_n));
        amp = fftshift( I(:, s(iterK)) );
        psi_n(Fsupp) = amp(Fsupp).*phase(Fsupp);
%         experimental soft thresholding, uncomment 2 following lines for try
%          psi_n(~Fsupp) = (real(psi_n(~Fsupp)) - del2 * sign(real(psi_n(~Fsupp)))).*(abs(psi_n(~Fsupp)) >= del2)+...
%                              1i*(imag(psi_n(~Fsupp)) - del2 * sign(imag(psi_n(~Fsupp)))).*(abs(psi_n(~Fsupp)) >= del2); 
        psi_n = ifft(psi_n)*N;
        
            Uo = conj(temp)./max( (abs(temp).^2) );
            Up = conj(Obj)./max( (abs(Obj).^2) );

            Corr1 = alpha.*Uo.*(psi_n - psi);
            Corr2 = sig_shift(alpha.*Up.*(psi_n - psi),-D(s(iterK)),F);

            Obj = Obj +  Corr1 + Corr2 ;
        Ir(:, s(iterK)) = abs( fftshift( fft(Obj.*temp)/N ) );
        
        
        
        
            if mod(iterK,K)== 0
                error(iter) = sqrt(sum(sum( abs(Ir(fftshift(Fsupp),:)-I(fftshift(Fsupp),:) ).^2 )))/sqrt(sum(sum( abs(I(fftshift(Fsupp),:) ).^2 )));
                fprintf('Iter:%d   IterK:%d alpha=%d Error=%d\n',iter, iterK, alpha, error(iter));
                
%                 subplot(2,2,1); 
%                 p1 = plot(time*1e15, abs(Obj), 'LineWidth',2);
%                 xlabel('Time [fsec]','FontSize',16); ylabel('Amplitude [a.u.]','FontSize',16);
% %                 xlim([-1e-13 1e-13]);
%                 title('Intensity Obj');
%                 
%                 subplot(2,2,2)
%                 p2 = plot(time*1e15, unwrap(angle(Obj))/pi, 'LineWidth',2);
%                 xlabel('Time [fsec]','FontSize',16); ylabel('Phase [Pi]','FontSize',16);
% %                 xlim([-1e-13 1e-13]);
%                 title('Phase Obj');
% 
%                 subplot(2,2,3)
%                 imagesc(time*1e15, fftshift(F)*1e-12, I.*kron(fftshift(Fsupp), ones(1, K)) );title('Used I');
%                 xlabel('Time [fsec]','FontSize',16); ylabel('Freq.[THz]','FontSize',16);
% 
%                 subplot(2,2,4)
%                 imagesc(time*1e15, fftshift(F)*1e-12, Ir);title('Recovered I');
%                 xlabel('Time [fsec]','FontSize',16); ylabel('Freq.[THz]','FontSize',16);
%                 pause(0.01);
            end
        
        
    end
    if error(iter)<STOPc
        return;
    end

    iter = iter+1;

end


    function [sig_out] = sig_shift(sig_in, d, F)
        sig_out = ifft( fft(sig_in).*exp(1i*2*pi*d*F) );
    end

end
