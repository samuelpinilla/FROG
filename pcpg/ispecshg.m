function [gpulse,ggate,error] = ispecshg(spectrogram,gpulse,ggate,A,iterations)

N = max(size(gpulse));
halfN = N/2;
temp = zeros(N,N);

for x=1:1:iterations
    efrog = gpulse'*ggate + conj(ggate')*conj(gpulse);
    
    for j=2:1:N
        temp(j,1:j-1) = efrog(j,1:j-1);
        efrog(j,1:N+1-j) = efrog(j,j:N);
        efrog(j,N+2-j:N) = temp(j,1:j-1);
    end
    
    temp(:,1:halfN) = efrog(:,1:halfN);
    efrog(:,1:halfN) = efrog(:,halfN+1:N);
    efrog(:,halfN+1:N) = temp(:,1:halfN);
    
    temp(1:halfN,:) = efrog(1:halfN,:);
    efrog(1:halfN,:) = efrog(halfN+1:N,:);
    efrog(halfN+1:N,:) = temp(1:halfN,:);
    
    efrog = fft(efrog);
    
    temp(1:halfN,:) = efrog(1:halfN,:);
    efrog(1:halfN,:) = efrog(halfN+1:N,:);
    efrog(halfN+1:N,:) = temp(1:halfN,:);
    
    for j=1:1:N
        for k=1:1:N
            temps = abs(efrog(j,k));
            if temps ~= 0
                efrog(j,k) = spectrogram(j,k)*(efrog(j,k)/temps);
            else
                efrog(j,k) = spectrogram(j,k);
            end
        end
    end
    
    temp(1:halfN,:) = efrog(1:halfN,:);
    efrog(1:halfN,:) = efrog(halfN+1:N,:);
    efrog(halfN+1:N,:) = temp(1:halfN,:);
    
    efrog = ifft(efrog);
    
    temp(1:halfN,:) = efrog(1:halfN,:);
    efrog(1:halfN,:) = efrog(halfN+1:N,:);
    efrog(halfN+1:N,:) = temp(1:halfN,:);
    
    temp(:,1:halfN) = efrog(:,1:halfN);
    efrog(:,1:halfN) = efrog(:,halfN+1:N);
    efrog(:,halfN+1:N) = temp(:,1:halfN);
    
    for j=2:1:N
        temp(j,N+2-j:N) = efrog(j,N+2-j:N);
        efrog(j,j:N) = efrog(j,1:N+1-j);
        efrog(j,1:j-1) = temp(j,N+2-j:N);
    end
    
    npulse = efrog'*gpulse';
    npulse = (efrog*npulse)';
    
    ggate  = efrog*ggate';
    ggate  = (efrog'*ggate)';
    
    gpulse = npulse/norm(npulse);
    ggate = ggate/norm(ggate);
    
    error(x)   = norm(spectrogram-abs(A(ggate.')),'fro')/norm(spectrogram,'fro');
    fprintf('Iter = %d  Error = %f \n',x, error(x));
    if min(error)<=1e-6
        break;
    end
end
return