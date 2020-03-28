function rho = coeff_circularity(v)
    c = mean(abs(v.^2));
    p = mean(v.^2);
    rho = abs(p)/c;
end