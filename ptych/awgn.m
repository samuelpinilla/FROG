function Ynoisy = awgn(in,snr,signalpower,seed)
    snraux = 10^(snr/10);
    sigma  = 1/snraux;
    rng(seed.Seed,seed.Type);
    noise  = in.*randn(size(in)).*sqrt(sigma);
    Ynoisy = in + noise;
end