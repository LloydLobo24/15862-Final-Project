function [pdK, pnIndex] = calculateCalibration(pdCalibrationLines, pdMask, nLineLength, nNumberCalibrationLines, nLeft, nRight)


pcdFFT = fft(pdCalibrationLines);

%pcdFFT = pcdFFT .* repmat(pdMask, [1 4]);
pcdFFT = pcdFFT.*pdMask;
pcdSpectrum = ifft(pcdFFT);

pdPhase = angle(pcdSpectrum);
pdPhase = unwrap(pdPhase, [], 1);

nMidPoint = (nLeft + nRight) / 2;
fTemp = pdPhase(nMidPoint, :);
pdPhase = pdPhase - repmat(fTemp, [nLineLength 1]);

dLeft = mean(pdPhase(nLeft, :));
dRight = mean(pdPhase(nRight, :));

dSlope = (dRight - dLeft) / (nRight - nLeft);
dOffset = dLeft - dSlope * nLeft;

pdK = zeros([nLineLength nNumberCalibrationLines]);
pnIndex = zeros([nLineLength nNumberCalibrationLines]);

for j = 1 : nNumberCalibrationLines
    
    pnAssigned = zeros([nLineLength, 1]);
    
    pdK(:, j) = (pdPhase(:, j) - dOffset) / dSlope;
    pnTemp = ceil(pdK(:, j));
    
    for i = 1 : nLineLength
        nTemp = pnTemp(i);
        if nTemp >= 1 && nTemp <= nLineLength
            pnIndex(nTemp, j) = i;
            pnAssigned(nTemp) = 1;
        end
    end
    
    i = 1;
    nTemp = 1;
    
    while pnAssigned(i) == 0 && i < nLineLength
        while i >= pdK(nTemp + 1, j) && nTemp < nLineLength
            nTemp = nTemp + 1;
        end
        
        pnIndex(i, j) = nTemp;
        i = i + 1;
    end
    
    if i <= nLineLength
        
        nTemp = i;
        nLast = pnIndex(i, j);
        if nLast < 1
            nLast = 1;
        end
        if nLast > nLineLength - 3
            nLast = nLineLength - 3;
        end
        pnIndex(i, j) = nLast;
        
        for i = nTemp + 1 : nLineLength
            if pnAssigned(i) == 1
                nLast = pnIndex(i, j);
                if nLast < 1
                    nLast = 1;
                end
                if nLast > nLineLength - 3
                    nLast = nLineLength - 3;
                end
                pnIndex(i, j) = nLast;
            else
                pnIndex(i, j) = nLast;
            end
        end
        
    else
        for i = 0 : nLineLength
            pnIndex(i, j) = i;
        end
    end
end

pnIndex = pnIndex - 1;

end