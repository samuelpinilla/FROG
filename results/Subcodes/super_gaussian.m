function out = super_gaussian(in , coeff_x, coeff_y, power)
   
    [n, m] = size(in);
    hx = -m/2: m/2-1;
    hy = -n/2: n/2-1;
    [X , Y] = meshgrid(hx,hy);
    wx = m/2 * coeff_x;
    wy = n/2 * coeff_y;
    sgauss = exp(-(X.^2/wx^2 + Y.^2/wy^2).^power);
    sgauss = sgauss + 1e-4; 
    sgauss = quickscale(sgauss);
    out = in .* sgauss;

end