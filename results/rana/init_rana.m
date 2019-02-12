function x1 = init_rana(Y,x)
%% initialization
N  = size(Y,1);
g  = sum(Y,2);

s = sqrt(ifft(ifftshift(g)));
s = s/norm(s,'fro');

%% variables
alpha = 0.09;
beta  = 0.425;
gamma = 1;

t0 = randsrc(1,1,3:N-3);
fprintf('t0=%d\n',t0);
ll = t0+1:1:N;

for r=1:length(ll)
    t = ll(r);
    
    term1p = alpha*abs(s(t)-s(t-1))^2;
    term1n = alpha*abs(-s(t)-s(t-1))^2;
    
    term2p = beta*abs(s(t)-s(t-1)-(s(t-1)-s(t-2)))^2;
    term2n = beta*abs(-s(t)-s(t-1)-(s(t-1)-s(t-2)))^2;
        
    term3p = gamma*abs(s(t)-s(t-1)-(s(t-1)-s(t-2)) - (s(t-1)-s(t-2)-(s(t-2)-s(t-3))))^2;
    term3n = gamma*abs(-s(t)-s(t-1)-(s(t-1)-s(t-2)) - (s(t-1)-s(t-2)-(s(t-2)-s(t-3))))^2;
    
    epsilon1 = term1p + term2p + term3p;
    epsilon2 = term1n + term2n + term3n;
    
    if epsilon1 > epsilon2
        s(t) = -s(t);
    end
end

ll = t0-1:-1:1;

for r=1:length(ll)
    t = ll(r);
    
    term1p = alpha*abs(s(t)-s(t+1))^2;
    term1n = alpha*abs(-s(t)-s(t+1))^2;
    
    term2p = beta*abs(s(t)-s(t+1)-(s(t+1)-s(t+2)))^2;
    term2n = beta*abs(-s(t)-s(t+1)-(s(t+1)-s(t+2)))^2;
        
    term3p = gamma*abs(s(t)-s(t+1)-(s(t+1)-s(t+2)) - (s(t+1)-s(t+2)-(s(t+2)-s(t+3))))^2;
    term3n = gamma*abs(-s(t)-s(t+1)-(s(t+1)-s(t+2)) - (s(t+1)-s(t+2)-(s(t+2)-s(t+3))))^2;
        
    epsilon1 = term1p + term2p + term3p;
    epsilon2 = term1n + term2n + term3n;
    
    if epsilon1 > epsilon2
        s(t) = -s(t);
    end
end

%% estimated spectrum

x0 = sqrt(fftshift(fft(s)/sqrt(N))).*exp(-1i*randsrc(N,1,[0,1]));
% x0 = x0/norm(x0);

%% real spectrum
% Re = fftshift(abs(fft(x)))/sqrt(N);

% error = norm(abs(x0)-Re)/N;
% fprintf('spectrum estimation: %f\n',error);

% figure;plot(1:N,abs(x0)/norm(x0),1:N,Re);

x1 = ifft(ifftshift(x0))*sqrt(N);
x1 = x1/norm(x1,'fro');
end