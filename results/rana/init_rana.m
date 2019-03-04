function x = init_rana(Isig,filter)
addpath(genpath('Subcodes'));
N  = size(Isig,1);
if filter
    Isig = super_gaussian(Isig,2,2,2);
    fI_shg = ifft_FROG(Isig);
    fI_shg = super_gaussian(fI_shg,100,.75,10);
    Isig = abs(fft_FROG(fI_shg));
    fI_shg = fftc(Isig, [], 2);
    fI_shg = super_gaussian(fI_shg,.75, 100,10);
    Isig = abs(ifftc(fI_shg,[],2));
    Isig  = remove_neg(Isig);
end

Isig = frogbacksub_noise(Isig, [5 5],[5 5]);
Mwin = sum(Isig,2);                                     %frequency marginal

Sout = SHG_marginal(Mwin, 1, 1) ;           % getting only one spectrum
x    = ifft(ifftshift(sqrt(Sout).*exp(1i*rand(1,N))));
x    = x.';
end
