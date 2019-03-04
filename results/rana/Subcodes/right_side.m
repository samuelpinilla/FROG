function [Swc, rms_err1, rms_err2] = right_side(Mw, St, pos_neg)


    N = length(pos_neg);
    
    signs = pos_neg;    
    signs(2:N/2) = fliplr(pos_neg(N/2+2:end));
 
    Stc = St.* signs;
    Swc = fftc(Stc);
    Swc = super_gaussian_1d(Swc, .65,5);

    A = abs(Swc);
    autoS = quickscale(abs(conv(A,  A,'same')));
 
    rms_err1 = mean(abs(autoS - Mw).^2); 
    % positivity
    [~ , ind]  = max(abs(real(Swc)));
    Swc = sign(real(Swc(ind))) * Swc; 
    rms_err2 = abs(sum(Swc(real(Swc) < 0)));
    Swc = abs(Swc);
 
end