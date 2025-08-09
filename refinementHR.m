function [isValidHR, estimateHR, prevValidHRs, countFrame] = refinementHR(segmentHRs, T1, T2, prevValidHRs, countFrame, Fs)

% Input parameters:
%   segmentHRs     — Heart rates of the three signal segments.
%   [T1, T2]             — High-probability interval for HR.
%   prevValidHRs    —  Heart rate sequence used for Case 2.
%
% Output parameters:
%   isValidHR           — Heart rate correction flag.
%   estimateHR        — the estimated heart rate.

    isValidHR = false; 
    estimateHR     = NaN;

    diffHR  = abs(diff(segmentHRs));
    meanHR    = mean(segmentHRs);

if ~any(diffHR > 22)
    if meanHR >= T1 && meanHR <= T2
        isValidHR    = true;
        estimateHR        = segmentHRs(end);

    elseif isempty(prevValidHRs)
        prevValidHRs        = segmentHRs;
        countFrame = -1;

    else
        diff2HR = abs(segmentHRs - prevValidHRs);
        if ~any(diff2HR > 10) && countFrame <= Fs
            isValidHR    = true;
            estimateHR        = segmentHRs(end);
        else
            prevValidHRs        = segmentHRs;
            countFrame = -1;
        end
    end
end

end
