function X = backsub_noise_1d(X, n1)
X = X(:)';
a = mean([X(1:n1) X(end-n1+1:end)]);
X = X - a;
X(X<0) = 1e-5;

end