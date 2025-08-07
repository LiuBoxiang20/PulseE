function [imfOptimized, alpha] = optimizeAlpha(data, Fs, freq, P1)

    ALPHA = 5000:10000:25000; 
    nA        = numel(ALPHA);
    tol        = 1e-7; 
    tau       = 0.001; 
    lenP1   = length(P1);

    ED   = zeros(1, nA);
    imfs = zeros(nA, numel(data));

    for k = 1:nA
        a = ALPHA(k);

        p = vme(data, a, freq, Fs, tau, tol);
        imfs(k, :) = p;
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
        ED(k) = sum((P1(ni)-P2(ni)).^2) + sum(P2(out_index).^2);
       
    end

    [~, lo] = min(ED);
    alpha = ALPHA(lo);
    imfOptimized  = imfs(lo, :);

end
