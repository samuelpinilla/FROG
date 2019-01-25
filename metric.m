function [val] = metric(x,z)
c1 = abs( xcorr(z,x) );
c2 = abs( xcorr(conj( z(end:-1:1) ),x) );
c3 = abs( xcorr(-z,x) );
c4 = abs( xcorr(conj( -z(end:-1:1) ),x)  );
[~, ind] = max([max(c1) max(c2) max(c3) max(c4)]);
switch ind
    case 1
        [~, indM] = max(abs(c1));
        pulseT = circshift( z, numel(z)-indM);
    case 2
        [~, indM] = max(abs(c2));
        pulseT = circshift( conj(z(end:-1:1)), numel(z)-indM);
    case 3
        [~, indM] = max(abs(c3));
        pulseT = circshift( -z, numel(z)-indM);
    case 4
        [~, indM] = max(abs(c4));
        pulseT = circshift( conj(-z(end:-1:1)), numel(z)-indM);
end
error = Inf;
N = numel(z);
F = ifftshift(-N/2:N/2-1);
for ii=1:400
    d = -1+ii/200;
    pulseR = ifft( fft(pulseT).*exp(1i*2*pi*d*F'/N) );
    if sqrt(sum(sum(abs( abs(pulseR)-abs(x) ).^2)))<error
        error = sqrt(sum(sum(abs( abs(pulseR)-abs(x) ).^2)));
        pulse = pulseR;
    end
end

%% metric
val = norm(x - exp(-1i*angle(trace(x'*pulse))) * pulse, 'fro')/norm(x,'fro');
end