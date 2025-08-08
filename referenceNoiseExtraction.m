function  refRawNoise = referenceNoiseExtraction(filePath, frameRange)

    gSignals = readmatrix(filePath);  % filePath contains the green-channel data for two background ROIs linked to the facial landmarks.
    gSignalsRange = gSignals(frameRange, :);  % frameRange matches the frame range used for extracting the rPPG signal.
    [~, score] = pca(gSignalsRange);
    refRawNoise =  score(:,1);

end

