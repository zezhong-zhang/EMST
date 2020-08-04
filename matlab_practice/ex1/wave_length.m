function lambda = wave_length(voltage)
    lambda = 12.2643/sqrt(voltage*10^3*(1+0.978476*10^-6*voltage*10^3));
end