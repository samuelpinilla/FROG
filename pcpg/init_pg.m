function g = init_pg(Y)
%% initialization
    N       = size(Y,1);

    g       = sum(Y.^2,1)'/sqrt(sum(sum(abs( sum(Y.^2,1) ).^2)));
    g       = g.*exp(-1i*randsrc(N,1,[0,1]));
end