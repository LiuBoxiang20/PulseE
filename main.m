clear; close all; clc; 

%% ==================== File Paths & Basic Parameters ======================
rppgDatafilePath             = '  ';  
referenceNoiseFilePath  =  '  ';  
exp_num                         = 1;            % Experiment number
subj_num                        = 1;            % Subject number
initialFrame                     = 1;            % Starting frame index
Fs                                    = 30;          % Sampling rate (frame rate) in Hz

convThr                           = 0.5;         % VPS convergence threshold
segmentIntervalSec        = 5;            % Time interval between segments (seconds)
l_l                                    = 0.75;       % Lower bound of heart rate band (Hz)
u_l                                   = 2.5;         % Upper bound of heart rate band (Hz)
countFrame                     = -1;
prevValidHRs                  = [];            % Store previously validated HR values
cycleIdx                           = 0;            % Iteration index
T1                                    = 64.6;       % [T1, T2] — High-probability interval for HR.
T2                                    = 93.6; 
NyquistF                          = Fs / 2; 
isValidHR                         = false;      % Heart rate correction flag.
SN                                    = Fs*20; 
l_f                                     = round(l_l*SN / Fs) + 1; 
u_f                                    = round(u_l*SN / Fs) + 1; 
bandIdx                            = l_f : u_f;  % Index range for HR frequency band
[B, A]                                = butter(4, [l_l/NyquistF, u_l/NyquistF]);  % 4th-order bandpass Butterworth filter

if isempty(gcp('nocreate'))
    parpool('local');  
end

%% ============================ Main Loop ==============================
while ~isValidHR
countFrame    = countFrame + 1;
cycleIdx          = cycleIdx + 1;
segmentHRs  = zeros(1, 3); 

parfor segmentIdx = 1:3
    
offset = (segmentIdx-1)*Fs*segmentIntervalSec + (cycleIdx-1);
segmentStart  = initialFrame + offset;

frameRange = segmentStart:segmentStart+Fs*20-1;  % 20-second segment

rppgRawData = rppgDataLoader(rppgDatafilePath, frameRange, Fs); rppgData = filtfilt(B, A, rppgRawData.').';  % Load rPPG data

refRawNoise = referenceNoiseExtraction(referenceNoiseFilePath, frameRange); refNoise = filtfilt(B,A,refRawNoise);  % Extract reference noise

Y = vps(rppgData, refNoise, Fs, convThr, bandIdx);  % rPPG signal decomposition and rigid motion noise removal

[rppgHR, selectedPulseSignal] = compositeIndex(Y, Fs, bandIdx); % Composite index: select the real pulse signal

segmentHRs(segmentIdx) = rppgHR;
end

 [isValidHR, estimateHR, prevValidHRs, countFrame] = refinementHR(segmentHRs, T1, T2, prevValidHRs, countFrame, Fs);  % Heart rate refinement

end

disp(['Estimated Heart Rate: ',num2str(estimateHR)])  % Output heart rate
