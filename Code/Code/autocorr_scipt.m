function [r_bias, r_unbias, lag] = autocorr_scipt(x)
    [r_bias, lag] = xcorr(x, 'bias');
    r_unbias = xcorr(x,'unbias');
end


