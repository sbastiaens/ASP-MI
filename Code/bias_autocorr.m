function [r_bias r_unbias] = bias_autocorr(x,N,K)
d=1;
for k=0:K
    r_tot = 0;
    for n = (k+1):N
       r = x(n)*conj(x(n-k));
       r_tot = r_tot+r;
    end
    r_bias(d,:) = (1/(N))*r_tot;
    r_unbias(d,:) = (1/(N-k))*r_tot;
    d=d+1;
end

end