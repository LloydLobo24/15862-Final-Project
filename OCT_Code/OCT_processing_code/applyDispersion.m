function pcdIMAQComplex = applyDispersion(pdIMAQ, pdR, pdI)

[~, nNumberLines] = size(pdIMAQ);

pdIMAQR = pdIMAQ .* repmat(pdR, [1 nNumberLines]);
pdIMAQI = pdIMAQ .* repmat(pdI, [1 nNumberLines]);
pcdIMAQComplex = complex(pdIMAQR, pdIMAQI);

clear pdIMAQR pdIMAQI

end