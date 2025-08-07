function Y = vps(rppgData, refNoise, Fs, convThr, bandIdx)

%     Input parameters:
%     rppgData  —  A raw (filtered) rPPG signal matrix of size [Nchan × Nsamp]
%     refNoise   —  A reference noise signal of length Nsamp
%     Fs            —  Sampling rate (Hz)  
%     convThr   —  VPS convergence threshold
%     maxIter    —  Maximum number of iterations per channel (e.g., 5)  
%
%     Output parameters:
%     Y              —  All extracted IMFs

    [Nchan, ~] = size(rppgData);
    Y = [];  
    maxIter = 5;

    for ch = 1:Nchan
        x0     = rppgData(ch, :);

        [~, ~, P1_0] = fftHR(x0, Fs);
        P0_S        = max(P1_0(bandIdx))^2;

        % iterative extraction
        residual = x0;
        imfs     = [];
        thr      = convThr;
        cnt      = 0;
        while cnt < maxIter
            cnt = cnt + 1;
            [~, freq, P1] = fftHR(residual, Fs);
            P_S        = max(P1(bandIdx))^2;
            if P_S / P0_S < thr
                break
            end

            [imfOptimized, alpha]  = optimizeAlpha(residual, Fs, freq, P1);
            residual = residual - imfOptimized;     
            imfOptimized     = rigidMotionRemoval(imfOptimized, refNoise, Fs, alpha);
            imfs     = [imfs; imfOptimized];                 
        end

        % If no IMF is extracted, further decompose the signal until a valid IMF is obtained.
        if isempty(imfs)
            while isempty(imfs)
                [~, freq, P1] = fftHR(residual, Fs);

                 [imfOptimized, alpha]  = optimizeAlpha(residual, Fs, freq, P1);
                residual = residual - imfOptimized;
                imfOptimized     = rigidMotionRemoval(imfOptimized, refNoise, Fs, alpha);
                imfs     = [imfs; imfOptimized];                
            end
        end

        Y = [Y; imfs];
    end
end
