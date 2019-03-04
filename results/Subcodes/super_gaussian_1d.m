function out = super_gaussian_1d(in , coeff, power)
    
    in = in(:)';
    [~, n] = size(in);
    h = -n/2: n/2-1;
       
    wx = n/2 * coeff;
    
    sgauss = exp(-(h.^2/wx^2 ).^power);
    
    sgauss = sgauss + 1e-4;
    sgauss = quickscale(sgauss);
   
    out = in .* sgauss;

end
