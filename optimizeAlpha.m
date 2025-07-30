function alpha = optimizeAlpha(data, Fs)

    ED = [];
    ALPHA = 5000:10000:25000; 
    tol = 1e-7; tau = 0.001; 
    [~, freq_s, P1] = fftHR(data, Fs); lenP1 = length(P1);

    for a = ALPHA
        p = vme(data, a, freq_s, Fs, tau, tol);
        [~, ~, P2] = fftHR(p, Fs); 
    
        [~, locs] = findpeaks(-P1);
        [~, peak_locs] = max(P1);
        distances = abs(locs-peak_locs);
        [~, sorted_indices] = sort(distances);
        nearest_indices = sort(locs(sorted_indices(1:2)));

        if peak_locs<nearest_indices(1) 
           [~, lo] = min(P1(1:peak_locs));
           nearest_indices = [lo, nearest_indices(1)];
        elseif peak_locs>nearest_indices(2)
            [~, lo] = min(P1(peak_locs:end));
           nearest_indices = [nearest_indices(2), lo + peak_locs - 1];
        end

        ni = nearest_indices(1):nearest_indices(2);
        out_index = [1:nearest_indices(1)-1 nearest_indices(2)+1:lenP1];
        ed = sum((P1(ni)-P2(ni)).^2) + sum(P2(out_index).^2);
        ED = [ED ed];

    end

    [~, lo] = min(ED);
    alpha = ALPHA(lo);

end
