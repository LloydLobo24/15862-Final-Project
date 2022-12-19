function pcdDepthProfile = getComplexDepthProfile(pcdIMAQ, pdMask)

[~, nNumberLines] = size(pcdIMAQ);

pcdIMAQMasked = pcdIMAQ .* repmat(pdMask, [1 nNumberLines]);
pcdDepthProfile = fft(pcdIMAQMasked);

end