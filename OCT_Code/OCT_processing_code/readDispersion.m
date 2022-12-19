function [pdR, pdI] = readDispersion(strFile)

fp = fopen(strFile, 'r', 'l');
nLineLength = fread(fp, 1, 'int');

pdR = fread(fp, nLineLength, 'float');
pdI = fread(fp, nLineLength, 'float');
fclose(fp);

clear fp nLineLength

end