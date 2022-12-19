clearvars; 
close all; 

%% parameters

% figure display parameters
scrsz = get(groot,'ScreenSize');
nY = 2;
nBottom = 50;
nTop = 90;
nX = 3;
nLeft = 10;
nRight = 10;
nHeight = scrsz(4)-nBottom;
nWidth = scrsz(3)-nLeft;
fA = figure('Position',[nLeft+0*nWidth/nX nBottom+nHeight/2 nWidth/nX-nRight nHeight/nY-nTop]);
fB = figure('Position',[nLeft+1*nWidth/nX nBottom+nHeight/2 nWidth/nX-nRight nHeight/nY-nTop]);
% fC = figure('Position',[nLeft+2*nWidth/nX nBottom+nHeight/2 nWidth/nX-nRight nHeight/nY-nTop]); 
% fD = figure('Position',[nLeft+0*nWidth/nX nBottom+0*nHeight/2 nWidth/nX-nRight nHeight/nY-nTop]);
% fE = figure('Position',[nLeft+1*nWidth/nX nBottom+0*nHeight/2 nWidth/nX-nRight nHeight/nY-nTop]);
% fF = figure('Position',[nLeft+2*nWidth/nX nBottom+0*nHeight/2 nWidth/nX-nRight nHeight/nY-nTop]);
clear scrsz nY nBottom nTop nX nLeft nRight nHeight nWidth;

% controling parameters
bMATLABInterp = true; 
bRef = 0;                % reference file: 1; no reference file: 0 
pnNumbers = 47; 
% nPos = pnNumbers(1); 

%% directories 
% data directories 
% strFolder = 'Sensitivity_20kHz\Sensitivity_Telecentric_1'; 
% strDir = ['D:\Raw data\USC System\20221118 Sensitivity\', strFolder, '\']; 
% strSaveMat = ['D:\Processed data\USC System\20220816 Phase stability\', strFolder, '\'];
% strSaveFig = ['D:\Processed data\USC System\20220816 Phase stability\', strFolder, '\figures\'];
% mkdir(strSaveMat); 
% mkdir(strSaveFig); 

strDir = 'G:\My Drive\Research\SD-OCT\120522\roll-off plot\'; 


listFolder = dir(sprintf('%s*.dat', strDir)); 



%% load and process common variables 
% load calibration and dispersion
strCalibration = 'G:\My Drive\Research\SD-OCT\OCT_Acquisition\nOCT 20210105\calibration_112922.dat';
% strCalibration = 'D:\Raw data\OCE and Scaffold\20220816 Demo\Calibration and dispersion\20220816_USC_calibration.dat';
[~, pdK, pnIndex] = readCalibration(strCalibration);
pdK = pdK(1, :)';
pnIndex = pnIndex(1, :)';

strDispersion = 'G:\My Drive\Research\SD-OCT\OCT_Acquisition\nOCT 20210105\dispersion_112922.dat'; 
% strDispersion = 'D:\Raw data\OCE and Scaffold\20220816 Demo\Calibration and dispersion\20220816_USC_dispersion.dat';
[pdDispReal, pdDispImag] = readDispersion(strDispersion); 


% load header from the first data file
strFile = fullfile(listFolder(1).folder, listFolder(1).name); 
try
    cellArrays = readHeader(strFile); 
catch
    cellArrays = readHeader2(strFile); 
end

nNumberLines = cellArrays{2,3};
nLineLength = cellArrays{2,4};

% calculate mask
nLeft = 256;
nRight = nLineLength - nLeft;
nRound = 32;
pdMask = calculateMask(nLineLength, nLeft, nRight, nRound);
 

% read reference spectrum / calculate reference spectrum
if bRef    % read reference spectrum directly from reference file      
    disp('reading reference...');
%     listRef = dir(sprintf('%s*Ref*.dat', strDir));
%     strReference = sprintf('%s%s', strDir, listRef(2).name);
    strReference = 'G:\My Drive\Research\SD-OCT\120522\roll-off plot\reference321435.dat';
    [pdIMAQ, ~] = readData(strReference, cellArrays); 
    pdRef = mean(pdIMAQ, 2);

    % calculate noise level
    pdRefSpectrum = pdIMAQ - repmat(pdRef, [1, size(pdIMAQ, 2)]); 
    clear pdIMAQ; 
    
    pdRefCalibrated = applyCalibration(pdRefSpectrum, pdK, pnIndex, bMATLABInterp);
    pcdRefCorrected = applyDispersion(pdRefCalibrated, pdDispReal, pdDispImag);
    pcdRefDepthProfiles = getComplexDepthProfile(pcdRefCorrected, pdMask);
    pcdRefDepthProfiles(round(nLineLength/2) + 1:end, : ) = []; 
    clear pdRefCalibrated pcdRefCorrected; 

    pdNoise = std(real(pcdRefDepthProfiles), 0, 2); 

    clear pcdRefDepthProfiles;

else    % calculate the reference spectrum by averaging spectra across multiple acquisition 
    disp('calculating reference...');
    pdRef = -1;

    % average multiple frames
%     for nNumber = pnNumbers(1)+ 20 : -1 : pnNumbers(1)+2
%         strFile = [listFolder(nNumber).folder, '\', listFolder(nNumber).name];
%         disp(strFile); 
% 
%         [pdIMAQ, ~] = readData(strFile, cellArrays); 
%         clear strFile; 
% 
%         if pdRef == -1
%             pdRef = mean(pdIMAQ, 2);
%         else
%             pdRef = 0.9*pdRef + 0.1*mean(pdIMAQ, 2);
%         end 
%         clear pdIMAQ;
%         
%     end
    clear nNumber; 

end

%% read actual data file
disp('processing data...');

% pddBSet = zeros(nLineLength/2, nNumberLines, length(pnNumbers)); 


for nNumber = pnNumbers

    % read IMAQ
    strFile = fullfile(listFolder(nNumber).folder, listFolder(nNumber).name); 
    fprintf('%d/%d: %s \n', nNumber-pnNumbers(1)+1, length(pnNumbers), strFile);
    [pdIMAQ, ~] = readData(strFile, cellArrays);

    %%% calculate reference spectrum based on single frame 
    pdRef = mean(pdIMAQ, 2);

    pdIMAQ = pdIMAQ - repmat(pdRef,[1, nNumberLines]);
    pdIMAQCalibrated = applyCalibration(pdIMAQ, pdK, pnIndex, bMATLABInterp);
    pcdIMAQCorrected = applyDispersion(pdIMAQCalibrated, pdDispReal, pdDispImag);
    pcdDepthProfiles = getComplexDepthProfile(pcdIMAQCorrected, pdMask);
    pcdDepthProfiles(round(nLineLength/2) + 1:end, : ) = []; 
    clear pdIMAQ pdIMAQCalibrated pcdIMAQCorrected; 

    % intensity image        
    pdI = abs(pcdDepthProfiles) .^ 2; 
    pddB = 10* log10(pdI); 

%     pdSNR = pdI ./ (repmat(pdNoise.^2, [1, size(pdI, 2)])); 


    % display 
    [~, strName, ~] = fileparts(strFile); 
    figure(fA), imagesc(pddB, [50, 110]), colormap(1-gray); colorbar;  
    title(sprintf('%d: %s', nNumber, strName), 'Interpreter','none'); 

    figure(fB), 
    plot(mean(pddB, 2));         
    grid on; xlim([1, 2048]); ylim([50, 110]); 
    xlabel('depth, pixel'); ylabel('intensity, dB');

    drawnow; 



end



