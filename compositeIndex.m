function [rppgHR, selectedPulseSignal] = compositeIndex(Y, Fs, bandIdx)

%   Input parameters:
%     Y    —  Candidate IMF matrix, with each row representing one IMF.
%     bandIdx —  Index range of the frequency interval.
%
%   Output parameters:
%     selectedPulseSignal   —  The selected pulse signal.
%     rppgHR                       —  The heart rate corresponding to the selected pulse signal.

    Ncomp = size(Y, 1);
    M = min(Ncomp, 6);
    indexAll = zeros(M, 5);

    % Compute each index and store the results in indexAll
    for i = 1:M
        candidate_imf = Y(i, :);
        [hrFFT, freq_s, P1_full] = fftHR(candidate_imf, Fs);
        P1 = P1_full(bandIdx);

        % Standard deviation of autocorrelation peak interval differences
        ac = xcorr(candidate_imf, 'biased');
        ac = ac / max(ac);
        [~, locs] = findpeaks(ac);
        addl = abs(diff(diff(locs))); 
        std_addl = std(addl);

        kurt_inv = 1 / kurtosis(P1);
        forfac = (impulsefactor(candidate_imf) / crestfactor(candidate_imf)) - 1.11;

        % Construct the index vector
        indexAll(i, :) = [std_addl, kurt_inv, forfac, hrFFT, freq_s];
    end

    index_n = size(indexAll, 2) - 2;    
    composite = sum(indexAll(:, 1:index_n), 2);
    normIndex = mapminmax(composite', 0, 1);
    [~, bestIdx] = min(normIndex);

    rppgHR = indexAll(bestIdx, end-1);
    selectedPulseSignal = Y(bestIdx, :);

end
