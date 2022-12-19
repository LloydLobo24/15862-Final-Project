function [cellArray, pdK, pnIndex] = readCalibration(strFile)

fp = fopen(strFile, 'r', 'l');
nNumberCalibrationLines = fread(fp, 1, 'int');
nNumberOCTLinesPerCalibration = fread(fp, 1, 'int');
nLineLength = fread(fp, 1, 'int');

cellArray(1) = nNumberCalibrationLines;
cellArray(2) = nNumberOCTLinesPerCalibration;
cellArray(3) = nLineLength;

pdK = fread(fp, nNumberCalibrationLines * nLineLength, 'float');
pdK = reshape(pdK, [nLineLength nNumberCalibrationLines])';

pnIndex = fread(fp, nNumberCalibrationLines * nLineLength, 'int');
pnIndex = reshape(pnIndex, [nLineLength nNumberCalibrationLines])';
fclose(fp);

clear fp nNumberCalibrationLines nNumberOCTLinesPerCalibration nLineLength

end