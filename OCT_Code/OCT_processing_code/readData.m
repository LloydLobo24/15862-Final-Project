function [pdIMAQ, pdDAQ] = readData(strFile, cellArrays)

fp = fopen(strFile, 'r', 'l');
for nArrayNumber = 2 : size(cellArrays, 1)
    strVar        = cellArrays{nArrayNumber, 1};
    nOffset       = cellArrays{nArrayNumber, 2};
    nNumberLines  = cellArrays{nArrayNumber, 3};
    nNumberPoints = cellArrays{nArrayNumber, 4};
    strDataType   = cellArrays{nArrayNumber, 5};
    fseek(fp, nOffset, 'bof');
    if strcmp(strVar, 'pdIMAQx2') % two cameras
        strTest = sprintf('%s = fread(fp, ''%s'');', strVar, strDataType);
        eval(strTest);        
        pdIMAQ1 = reshape(pdIMAQx2(1:2:end), [nNumberPoints, nNumberLines]); % pdIMAQParallel (H)
        pdIMAQ2 = reshape(pdIMAQx2(2:2:end), [nNumberPoints, nNumberLines]); % pdIMAQPerpendicular (V)
        pdIMAQ(:,:,1) = pdIMAQ1; 
        pdIMAQ(:,:,2) = pdIMAQ2;        
        clear pdIMAQ1 pdIMAQ2; 
    else
        strTest = sprintf('%s = fread(fp, [%d, %d], ''%s'');', strVar, nNumberPoints, nNumberLines, strDataType);
        eval(strTest);
    end    
    clear strVar nOffset nNumberLines nNumberPoints strDataType;
    clear strTest;
end
fclose(fp);
clear ans fp nArrayNumber;

end