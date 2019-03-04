function [Swc, rms_err1, rms_err2] = left_side(Mw, St, pos_neg)


    N = length(pos_neg);
    
    signs = pos_neg;    
    signs(N/2+2:end) = fliplr(pos_neg(2:N/2));

 
    Stc = St.* signs;
    Swc = fftc(Stc);
    Swc = super_gaussian_1d(Swc, .65,5);
    
    A = abs(Swc);  
    autoS = quickscale(abs(conv(A, A,'same'))); 
    rms_err1 = mean(abs(autoS- Mw).^2); 
    
    % positivit
    [~ , ind]  = max(abs(real(Swc)));
    Swc = sign(real(Swc(ind))) * Swc; 
    rms_err2 = abs(sum(Swc(real(Swc) < 0)));
    Swc = abs(Swc);
 
end