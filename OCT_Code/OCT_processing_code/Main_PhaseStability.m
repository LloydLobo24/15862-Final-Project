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
fC = figure('Position',[nLeft+2*nWidth/nX nBottom+nHeight/2 nWidth/nX-nRight nHeight/nY-nTop]); 
% fD = figure('Position',[nLeft+0*nWidth/nX nBottom+0*nHeight/2 nWidth/nX-nRight nHeight/nY-nTop]);
fE = figure('Position',[nLeft+1*nWidth/nX nBottom+0*nHeight/2 nWidth/nX-nRight nHeight/nY-nTop]);
% fF = figure('Position',[nLeft+2*nWidth/nX nBottom+0*nHeight/2 nWidth/nX-nRight nHeight/nY-nTop]);
clear scrsz nY nBottom nTop nX nLeft nRight nHeight nWidth;

% controling parameters
bMATLABInterp = true; 
bRef = 1;                % reference file: 1; no reference file: 0 
pnNumbers = 5:12; 
pnSets = 2 : 9; 


% physical parameters 
dLineRate = 20e3;   % Hz


%% directories 
% data directories 
strFolder = 'PhaseStability_20kHz'; 
strDir = ['D:\Raw data\USC System\20220816 Phase stability\', strFolder, '\']; 
strSaveMat = ['D:\Processed data\USC System\20220816 Phase stability\', strFolder, '\'];
% strSaveFig = ['D:\Processed data\USC System\20220816 Phase stability\', strFolder, '\figures\'];
mkdir(strSaveMat); 
% mkdir(strSaveFig); 

listFolder = dir(sprintf('%s*.dat', strDir)); 

%% load and process common variables 
% load calibration and dispersion
strCalibration = 'D:\Calibration and Dispersion files\USC\20220816_USC_calibration_2.dat';
% strCalibration = 'D:\Raw data\OCE and Scaffold\20220816 Demo\Calibration and dispersion\20220816_USC_calibration.dat';
[~, pdK, pnIndex] = readCalibration(strCalibration);
pdK = pdK(1, :)';
pnIndex = pnIndex(1, :)';

strDispersion = 'D:\Calibration and Dispersion files\USC\20220816_USC_dispersion_2.dat'; 
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

pdCompositeSNR = zeros(length(pnSets), 1); 
pdSigma = zeros(length(pnSets), 1);
nCounter = 1; 

for nSet = pnSets

    % read reference spectrum / calculate reference spectrum
    if bRef         
        disp('reading reference...');
        listRef = dir(sprintf('%s*NDF_%d_Ref*.dat', strDir, nSet));
        strReference = sprintf('%s%s', strDir, listRef(2).name);
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

%         pdNoise = mean(abs(pcdRefDepthProfiles).^2 , 2); 
%         pdNoise2 = std(pcdRefDepthProfiles, 0, 2); 
        pdNoise = std(real(pcdRefDepthProfiles), 0, 2); 
%         pdNoise4 = std(imag(pcdRefDepthProfiles), 0, 2);

%         figure, p1 = plot(10*log10(pdNoise)); 
%         hold on; p2 = plot(20*log10(pdNoise2)); hold off;  
%         hold on; p3 = plot(20*log10(pdNoise3)); hold off;  
%         hold on; p4 = plot(20*log10(pdNoise4)); hold off;         
%         grid on; xlim([1, 2048]); 
%         legend([p1 p2], {'mean(abs(pdComplex)^2), old method', 'std(real(pdComplex))'}); 

    % 
    %         clear pcdRefDepthProfiles;
    
    else 
        disp('calculating reference...');
        pdRef = -1;
        for nNumber = pnNumbers(1)+ 20 : -1 : pnNumbers(1)+2
            strFile = [listFolder(nNumber).folder, '\', listFolder(nNumber).name];
            disp(strFile); 
    
            [pdIMAQ, ~] = readData(strFile, cellArrays); 
            clear strFile; 
    
            if pdRef == -1
                pdRef = mean(pdIMAQ, 2);
            else
                pdRef = 0.9*pdRef + 0.1*mean(pdIMAQ, 2);
            end 
            clear pdIMAQ;
            
        end
        clear nNumber; 
    
    end
    
    %% read actual data file
    disp('processing data...');
    
    pddBSet = zeros(nLineLength/2, nNumberLines, length(pnNumbers));     
    
    pcdComplex = []; 
    for nNumber = pnNumbers
         
        % read IMAQ 
        listSet = dir(sprintf('%s*NDF_%d*.dat', strDir, nSet));
        strFile = fullfile(listSet(nNumber).folder, listSet(nNumber).name); 
        fprintf('%d/%d: %s \n', nNumber-pnNumbers(1)+1, length(pnNumbers), strFile);
        [pdIMAQ, ~] = readData(strFile, cellArrays);
    
        pdIMAQ = pdIMAQ - repmat(pdRef,[1, nNumberLines]);
        pdIMAQCalibrated = applyCalibration(pdIMAQ, pdK, pnIndex, bMATLABInterp);
        pcdIMAQCorrected = applyDispersion(pdIMAQCalibrated, pdDispReal, pdDispImag);
        pcdDepthProfiles = getComplexDepthProfile(pcdIMAQCorrected, pdMask);
        pcdDepthProfiles(round(nLineLength/2) + 1:end, : ) = []; 
        clear pdIMAQ pdIMAQCalibrated pcdIMAQCorrected; 

        pcdComplex = [pcdComplex, pcdDepthProfiles]; 

%         keyboard; 
        [~, strName, ~] = fileparts(strFile); 
        save(sprintf('%s%s.mat', strSaveMat, strName), 'pcdDepthProfiles', 'pdNoise'); 

    end

    % actual time 
    pdX = (1 : size(pcdComplex, 2))'; 
    pdT = (pdX - 1) / dLineRate; 

    % intensity image        
    pdI = abs(pcdComplex) .^ 2; 
    pddB = 10* log10(pdI); 

%         % processing: surface
%         pddBFilt = imfilter(pddB, ones([4, 16])/(4*16), 'replicate');
%         pnSurface = zeros(1, size(pdI, 2)); 
%     
%         [pnX, pnY] = find(pddBFilt(51:end, :) > 80);
%         pnSurface(flipud(pnY)) = flipud(pnX) + 50;
%         clear pnX pnY;
%     
%         % align with surface 
%         nSurf = round(mean(pnSurface(101:400))); 
%         nShift = -(nSurf - 200); 
%         pddB2 = circshift(pddB, [nShift, 0]); 
%     
%         pddBSet(:, :, nCounter) = pddB2; 
%     
%         clear nSurf nShift pddB2; 

    % processing: phase vs SNR
    pdSNR = pdI ./ (repmat(pdNoise.^2, [1, size(pdI, 2)])); 

    % phase difference
    pdSNRLine = mean(10*log10(pdSNR), 2); 
    [~, n1] = max(pdSNRLine(180:240)); n1 = n1 + 180 - 1; 
    [~, n2] = max(pdSNRLine(340:400)); n2 = n2 + 340 - 1; 
    nPhaseLine = [n1, n2]; 

    ddBNoise2 = mean(pdSNRLine(n1+20 : n2-20)); 
    pddBNoise2(nCounter) = ddBNoise2; 


    
    pdPhaseDiff = processPhaseDiff(pcdComplex, nPhaseLine);
    pdPhaseDiff = highpass(pdPhaseDiff, 200, dLineRate); 
%     pdPhaseDiff = bandstop(pdPhaseDiff, [50,70], dLineRate); 
    pdPhaseFFT = fftshift(fft(pdPhaseDiff));
    dSigma = std(pdPhaseDiff); 
    pdSigma(nCounter) = dSigma; 

    % SNR
    dSNR1 = pdSNR(n1); 
    dSNR2 = pdSNR(n2); 
    dSNR = 1 / (0.5/dSNR1 + 0.5/dSNR2); 
    pdCompositeSNR(nCounter) = dSNR; 

    % phase difference fft analysis
     
    dF = dLineRate / length(pdPhaseDiff); 
    pdFreq = (-dLineRate/2 : dF : dLineRate /2 - dF)'; 
    clear dF; 
    

    % display 
    [~, strName, ~] = fileparts(strFile); 
    figure(fA), imagesc(pddB(1:1024, :), [50, 100]), colormap(1-gray); colorbar; 
%         hold on; plot(pnSurface, 'y', 'LineWidth', 2); hold off; 
    title(strName, 'Interpreter','none'); 

    figure(fB), 
    subplot(2, 1, 1); 
    p1 = plot(mean(pddB, 2)); 
    hold on; p2 = plot(20*log10(pdNoise)); hold off;         
    grid on; xlim([1, 2048]); ylim([55, 120]); 
    legend([p1 p2], {'coverslip data', 'reference noise, std(real(pdComplex))'});
    title('coverslip and reference profiles')
    subplot(2, 1, 2); 
    p3 = plot(pdSNRLine); 
    hold on; plot(n1, pdSNRLine(n1), 'r*'); 
    hold on; plot(n2, pdSNRLine(n2), 'r*'); hold off;  
    grid on; xlim([1, 2048]); ylim([0, 55])
    title('SNR, (|E|/|n|)^2'); 

    figure(fC); 
    subplot(2, 1, 1); 
    plot(pdT*1000, pdPhaseDiff);     
    xlim([0, max(pdT*1000)]); ylim([-4, 4]); 
    grid on; 
    xlabel('Time, ms'); ylabel('\Delta\phi, rad'); 
    title('phase difference')

    subplot(2, 1, 2); 
    plot(pdFreq, 10*log10(abs(pdPhaseFFT))); 
    xlim([0, 1000]); ylim([-10, 25]);
    grid on; 
    xlabel('Frequency, Hz'); 
    title('FFT of phase difference')

    drawnow; 

    nCounter = nCounter + 1; 
%     keyboard; 


end 

% convert scales 
pddBSNR = 10 * log10(pdCompositeSNR) - pddBNoise2'; 

x = 1 : 50; 
xSNR = 10.^(x/10); 
yExpectedLine = sqrt(1./xSNR); 
figure(fE), 
p1 = plot(pddBSNR, pdSigma, '*r'); 
hold on; p2 = plot(x, yExpectedLine, 'b'); 
hold off; hold off; hold off; 
xlim([5, 50]); ylim([0.001, 5]); xlabel('SNR (dB)'); ylabel('standard deviation'); 
legend([p1, p2], {'measured', 'theoretical'}, 'Location', 'northeast'); 
set(gca, 'YScale', 'log');  grid on; 
title('phase stability'); 

%
% pdExpectedSNR = 10.^(pddBSNR/10);
% pdExpectedSigma = sqrt(1./pdExpectedSNR); 
% figure, 
% p1 = plot(pddBSNR, pdSigma, '*r'); 
% hold on; p2 = plot(pddBSNR, pdExpectedSigma, 'ob'); hold off; 
% hold on; p3 = plot(pddBSNR, pdSigma - pdExpectedSigma); hold off; 
% xlim([5, 50]); ylim([0, 0.5]) xlabel('SNR (dB)'); ylabel('standard deviation'); 
% legend([p1, p2, p3], {'measured', 'theoretical', 'measured - theoretical'}, 'Location', 'northeast'); 
% grid on; 




function pdPhaseDiff = processPhaseDiff(pdDepthProfiles, nPhaseLines)

nTop    = nPhaseLines(1); 
nBottom = nPhaseLines(2); 

pdAngleTop      = unwrap(angle(pdDepthProfiles(nTop, :))); 
pdAngleBottom   = unwrap(angle(pdDepthProfiles(nBottom, :))); 
pdRawPhaseDiff  = pdAngleTop - pdAngleBottom; 
if abs(mean(pdRawPhaseDiff)) > 3.0
    pdRawPhaseDiff  = pdRawPhaseDiff + pi/4;
    pdPhaseDiff     = wrapToPi(pdRawPhaseDiff) - pi/4; 
else
    pdPhaseDiff     = wrapToPi(pdRawPhaseDiff); 
end
clear nTop nBottom pdAngleTop pdAngleBottom; 

end

