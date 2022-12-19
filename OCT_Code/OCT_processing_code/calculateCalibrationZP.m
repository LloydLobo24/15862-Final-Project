  function [pdCalibrationZPPoint,pnCalibrationIndex] = calculateCalibrationZP(pdCalibrationLines, pdMask, nLineLength,nNumberLines,nSkipLines,nZPFactor,)

  nZPLineLength = nLineLength*nZPFactor;
  
 
  
  if(nSkipLines<1)
      nSkipLines = 1;
  end
  if(nSkipLines>64)
      nSkipLines = 64;
  end
  
  
  dLineCount = 0;
  for nLine = 1:nSkipLines:nNumberLines
      
      pcdFFT = fft(pdCalibrationLines);
      pdSum = 
      pcdMaskedFFT = pcdFFT.*pdMask;
      
      pcdSpectrum = ifft(pcdMaskedFFT);
      pdPhase = angle(pcdSpectrum);
      pdPhase = unwrap(pdPhase, [], 1);
      
      pdZPSum = 
      dLineCount = dLineCount+1;
      
  end
  
  dMin = pdZPSum(1)/dLineCount;
  dMax = pdZPSum(nZPLineLength)/dLineCount;
  
  dSlope = (dMax-dMin)/nZPLineLength;
  
  for nPoint = 1:nZPLineLength
  
      pdCalibrationPhase(1,nPoint) = pdZpSum(nPoint)/dLineCount;
      pdCalibrationPhase(2,nPoint) = dMin + dSlope * nPoint;
      pdCalibrationZPPoint(nPoint) = nPoint + (pdCalibrationPhase(0,nPoint)-pdCalibrationPhase(1,nPoint))/dSlope;
      
      
  end
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
end

