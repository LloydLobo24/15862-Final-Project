function cellArrays = readHeader2(strFile)

fp = fopen(strFile, 'r', 'n', 'US-ASCII');

% read filename
%nLength = fread(fp, 1, 'int');  strTest = transpose(char(fread(fp, nLength+1, 'char')));  strTest(1) = [];
nLength = fread(fp, 1, 'int');  strTest = transpose(char(fread(fp, nLength+2, 'char')));  strTest(1:2) = [];
strFilename = strTest;
clear nLength strTest;
% read frame number

nLength = fread(fp, 1, 'int');  strTest = transpose(char(fread(fp, nLength+1, 'char')));  strTest(1) = [];  eval(strTest);
clear nLength strTest;
% read number of data arrays

nLength = fread(fp, 1, 'int');  strTest = transpose(char(fread(fp, nLength+1, 'char')));  strTest(1) = [];  eval(strTest);
clear nLength strTest;
cellArrays{nNumberDataArrays+1, 5} = 0;
cellArrays{1,1} = strFilename;
cellArrays{1,2} = nFrameNumber;
clear strFilename nFrameNumber;
for nArrayNumber = 1 : nNumberDataArrays
    % read variable name, offset, number of lines, number of points, data type
    nLength = fread(fp, 1, 'int');  strTest = transpose(char(fread(fp, nLength+1, 'char')));  strTest(1) = [];  eval(strTest);  clear nLength strTest;
    nLength = fread(fp, 1, 'int');  strTest = transpose(char(fread(fp, nLength+1, 'char')));  strTest(1) = [];  eval(strTest);  clear nLength strTest;
    nLength = fread(fp, 1, 'int');  strTest = transpose(char(fread(fp, nLength+1, 'char')));  strTest(1) = [];  eval(strTest);  clear nLength strTest;
    nLength = fread(fp, 1, 'int');  strTest = transpose(char(fread(fp, nLength+1, 'char')));  strTest(1) = [];  eval(strTest);  clear nLength strTest;
    nLength = fread(fp, 1, 'int');  strTest = transpose(char(fread(fp, nLength+1, 'char')));  strTest(1) = [];  eval(strTest);  clear nLength strTest; 
    cellArrays{nArrayNumber+1, 1} = strVar;
    cellArrays{nArrayNumber+1, 2} = nOffset;
    cellArrays{nArrayNumber+1, 3} = nNumberLines;   % previously nNumberPoints
    cellArrays{nArrayNumber+1, 4} = nLineLength;
    cellArrays{nArrayNumber+1, 5} = strDataType;
    clear nOffsetIMAQ nNumberLines nNumberPoints strDataType;
end % for nArrayNumber

clear nArrayNumber nNumberDataArrays;
fclose(fp);
clear fp ans; 

end
