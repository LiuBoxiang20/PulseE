function IMFs_removed = rigidMotionRemoval(IMFs, ref_noise, Fs, alpha)

    tol = 1e-7; tau = 0.001; 
    EID = false(1, size(IMFs,1));

    for i = 1:size(IMFs, 1)
        [~, freq] = fftHR(IMFs(i,:), Fs);

        noiseIMF = vme(ref_noise, alpha, freq, Fs, tau, tol);
        [correlation, lag] = xcorr(noiseIMF, IMFs(i, :));
        [~, max_index] = max(abs(correlation));
        delay = lag(max_index);

        EID(i) = abs(delay) <= 0.1 * Fs;
    end

        IMFs_removed = IMFs(~EID, :);    

end
    