function out = remove_neg(in)

out = in ;

out(in<0) = 1e-50;

end

