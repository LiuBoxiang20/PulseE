function rppgRawData = rppgDataLoader(filePath, frameRange, Fs)

%{
Inputs:
  filePath          — Path to a matrix file containing the raw rPPG signals obtained by 
                             performing frame-wise spatial averaging over each facial ROI. 
                             If there are 3 ROIs, columns 1–3 hold the R/G/B channels for ROI 1, 
                             columns 4–6 for ROI 2, and columns 7–9 for ROI 3, and so on.
  frameRange  — Vector of frame indices to retain (row numbers).
  Fs                  — Frame rate in Hz.

Output:
  rppgRawData — 2×M matrix, where
                              row 1 is the raw rPPG signal extracted by CHROM,
                              row 2 is the raw rPPG signal extracted by POS.
%}

    rppgSignals = readmatrix(filePath); 
    rppgSignals = rppgSignals(frameRange, :);
   
    [~, C] = size(rppgSignals);
    rCols = 1:3:C;
    gCols = 2:3:C;
    bCols = 3:3:C;
    meanR = mean(rppgSignals(:, rCols), 2);   
    meanG = mean(rppgSignals(:, gCols), 2);  
    meanB = mean(rppgSignals(:, bCols), 2);   
    RGB = [meanR meanG meanB];
    
    N = size(RGB,1);                
    bvpRawCHROM = zeros(1,N);    
    bvpRawPOS = zeros(1,N);    
    windowDuration = 1.2;
    wl = ceil(windowDuration*Fs);      

    for start = 1:(N - wl + 1)
        idx = start : (start + wl - 1);       
        Cn = ( RGB(idx,:) ./ mean(RGB(idx,:)) )'; 

        % -------- CHROM --------
        S1  = [3, -2, 0; 1.5, 1, -1.5] * Cn;               
        bvpWin1  = S1(1,:) - (std(S1(1,:))/std(S1(2,:))) * S1(2,:);  
        bvpWin1  = bvpWin1 - mean(bvpWin1); 
        bvpRawCHROM(idx) = bvpRawCHROM(idx) + bvpWin1;       
        
        % -------- POS --------
        S2  = [0, 1, -1; -2, 1, 1] * Cn;               
        bvpWin2  = S2(1,:) + (std(S2(1,:))/std(S2(2,:))) * S2(2,:);  
        bvpWin2  = bvpWin2 - mean(bvpWin2);       
        bvpRawPOS(idx) = bvpRawPOS(idx) + bvpWin2;                
    end

    rppgRawData = [ bvpRawCHROM;  bvpRawPOS ];

end

