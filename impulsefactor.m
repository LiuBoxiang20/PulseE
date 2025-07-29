function pulsef = impulsefactor(x)
    pulsef = peak(x)/mean(abs(x));
end