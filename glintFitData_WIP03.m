clear all
dbstop if error

addpath('/Users/bnorris/code/GLINT/fitChiSquare') %For the fitchisq code
addpath('fitChiSquare') %For the fitchisq code

%dataDir = './';
% dataDir = '/Users/bnorris/Don''t Backup/NullerDataTemp/BinnedSaved/';
%dataFileName = 'BinnedData_alfBooGoodfiles_WithErrs_20160319T032628-20160319T034741_binsize100_peakEstBinSize100';
%dataFileName = 'BinnedData_Labtest200nm_WithErrs_20150323T001134-20170323T001648_binsize100_peakEstBinSize100';
%dataFileName = 'BinnedData_alfHer_WithErrs_20160321T041027-20160321T044159_binsize100_peakEstBinSize100';
%dataFileName = 'BinnedData_alfLyn_WithErrs_20160318T235416-20160319T000216_binsize100_peakEstBinSize100';
%dataFileName = 'BinnedData_Altair_WithErrs_20160320T055046-20160320T061028_binsize100_peakEstBinSize100';

% dataDir = '/Volumes/silo4/snert/GLINT_Data/2016AugustSubaru/BinnedSaved/';
% dataFileName = 'BinnedData_20160814T224045-20160814T234820_binsize100_peakEstBinSize100'

% dataDir = '/Users/bnorris/DontBackup/GLINTdata/20170531_Subaru/BinnedSaved/';
%dataDir='';
dataDir='..\GLINT_TestData\';
dataFileName = 'BinnedData_20170531T051354-20170531T052604_binsize100_peakEstBinSize100_NSC';
dataFileName = 'BinnedData_20160319T032628-20160319T034741_binsize100_peakEstBinSize100_NSC';
% dataFileName = 'BinnedData_20160319T032628-20160319T034741_binsize100_NSC';

dataDir='..\GLINT_BinnedSaved\';
% dataDir='..\GLINT_BinnedSaved\2017May-BinnedSaved\';
% dataFileName = 'BinnedData_20160319T032628-20160319T034741_binsize100_peakEstBinSize100_NSC_maxHist5'
% dataFileName = 'BinnedData_20160319T032628-20160319T034741_binsize100_peakEstBinSize100_NSC'
% dataFileName = 'BinnedData_20160319T032628-20160319T034741_binsize100_NSC'
dataFileName = 'BinnedData_alfBoo-SubMar2016_fixedPhotGT0_20160319T032628-20160319T034741_binsize100_NSC_peakest2';
dataFileName = 'BinnedData_alfBoo-SubMar2016_20160319T032628-20160319T034741_binsize100_peakEstBinSize100_NSC';

dataFileName = 'BinnedData_betPeg-SubMay2017-bestfewfiles_20170531T051354-20170531T052604_binsize100_peakEstBinSize100_NSC';
% dataFileName = 'BinnedData_betPeg-SubMay2017-bestfewfiles_20170531T051354-20170531T052604_binsize100_NSC_peakest2'

dataFileName = 'BinnedData_20160815T021820-20160815T022026_binsize100_peakEstBinSize100_NSC';

fitSaveFileName = 'fittedParams_wip_test';

% dataFileName = 'BinnedData_omiCet_subset1_SubNov2016_20161109T022345-20161109T023622_binsize100_peakEstBinSize100_NSC'
% fitSaveFileName = 'fittedParams_wip99';
%fitSaveFileName = 'fittedParams_20170804a';


global nSamps nLoops histEdgesSpecify guessParams gpuMode NSCPDFs


% Basin hopping
numHops = 1;%100;
hopSigmas = [0.1 1 0.01 0.1 1 0.01 0.01 0.01 0.01 0.01];
% hopSigmas = hopSigmas*10
hopSigmas = hopSigmas/2
% hopSigmas = hopSigmas/5
% hopSigmas = zeros(1,10)
relStepSigma = 1; %finDiffRelStep scales as 10^N(0,relStepSigma);


fitAll = false; % false to just fit the 5 free params
performFit = true;
fitUncertainties = false;
ignoreErrors = false;
figNum = 1

histPlotYlim = [0, 5]; %For plotting

% Use NSC? If false, use ASC.
useNSC = true;

gpuMode = true; % Set to true to use GPU based MC function (WIP)

% MC options
nLoops=256;%16;%256
nSamps = 2^16;

nLoops=128;
nSamps = 2^24;

nLoops=1;%8;%16;%256
nSamps = 2^16;%2^24;

nLoops=32;%128;
nSamps = 2^24;

% Restrict fitting to a smaller domain (e.g. to exclude long tail of zeros)
restrictedFittingDomain = [];
restrictedFittingDomain = [-1., 1.2];
restrictedFittingDomain = [-0.2, 1.2];
% restrictedFittingDomain = [-0.5, 2];

% Set up pdf measurement options
estimateHistErrors = true; % If true, use the Bernoulli RV method
% histNBins = 100;
% minXVal = -0.2;
% maxXVal = 1.2;
histNBins = 200;
minXVal = -1.2;
maxXVal = 1.2;
% minXVal = -1.5;
% maxXVal = 1.5;

finDiffRelStep = 1e-2; %1e-2 seems to work well

% Error scaling:
% Give priority to low-value measurements, as unmodelled errors (eg
% correlation between phi and I) are worse at high values.
errorScalingCutoff = 0.2;
errorScalingFactor = 1;


% Specify Guess parameters
deltaPhiMu = 0.3;
deltaPhiSig = 0.5*pi;
astroNull = 0.02;

% Ir would ideally not be free.
% Ir is the relative intensity deviation. In Hanot, etc, this is largely to
% account for non-simultaneous measurement of intensity. In our case, maybe we
% care more that I+ is a bad measurement, since assumes dPhi=0... ?
IrMu = 1.;
IrSig = 0.01;


% Fiddling for Subaru May2017 data
% deltaPhiMu = 0.3;
% deltaPhiSig = 0.4*pi;
% astroNull = 0.05;
% IrMu = 1.;
% IrSig = 0.8;
deltaPhiMu = 0.3;
deltaPhiSig = 0.3*pi;
astroNull = 0.02;
IrMu = 1.;
IrSig = 0.01;
% astroNull = 0.05;


% % Starting near (2 sig fig) best params from multi-hop fit
% deltaPhiMu = -1.566
% deltaPhiSig = -0.11006
% astroNull = -0.445
% IrMu = 0.96859
% IrSig = 0.0018867


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totsamps = nSamps * nLoops;
disp(['Total samps: ' num2str(totsamps, '%.2g')])

if useNSC
    load([dataDir dataFileName], 'nullEst', 'nullEstErr', 'measuredStats', ...
        'binCents', 'nullEstPDFVals', 'phot1DarkPDFVals', 'phot2DarkPDFVals', ...
        'nullChanDarkPDFVals', 'brightChanDarkPDFVals', 'deltaIPDFVals', ...
        'totalIPDFVals', 'phot1PDFVals', 'phot2PDFVals');
    NSCPDFs = struct('binCents', binCents, 'nullEstPDFVals', nullEstPDFVals, ...
        'phot1DarkPDFVals', phot1DarkPDFVals, 'phot2DarkPDFVals', phot2DarkPDFVals, ...
        'nullChanDarkPDFVals', nullChanDarkPDFVals, ...
        'brightChanDarkPDFVals', brightChanDarkPDFVals, ...
        'deltaIPDFVals', deltaIPDFVals, 'totalIPDFVals', totalIPDFVals, ...
        'phot1PDFVals', phot1PDFVals, 'phot2PDFVals', phot2PDFVals);
else
    load([dataDir dataFileName], 'nullEst', 'nullEstErr', 'measuredStats', ...
        'binCents', 'nullEstPDFVals');
    NSCPDFs = [];
end
nullEstBinCents = binCents;
nullEstPDFValsErr = zeros(1,length(nullEstPDFVals));

measuredNullEst = nullEst;
measuredNullEstErr = nullEstErr;
clear nullEst nullEstErr


if ~isempty(restrictedFittingDomain)
    if useNSC
        NSCInds = (NSCPDFs.binCents > restrictedFittingDomain(1)) & ...
            (NSCPDFs.binCents < restrictedFittingDomain(2));
        NSCPDFs.binCents = NSCPDFs.binCents(NSCInds);
        NSCPDFs.nullEstPDFVals = NSCPDFs.nullEstPDFVals(NSCInds);
        NSCPDFs.phot1DarkPDFVals = NSCPDFs.phot1DarkPDFVals(NSCInds);
        NSCPDFs.phot2DarkPDFVals = NSCPDFs.phot2DarkPDFVals(NSCInds);
        NSCPDFs.nullChanDarkPDFVals = NSCPDFs.nullChanDarkPDFVals(NSCInds);
        NSCPDFs.brightChanDarkPDFVals = NSCPDFs.brightChanDarkPDFVals(NSCInds);
        NSCPDFs.deltaIPDFVals = NSCPDFs.deltaIPDFVals(NSCInds);
        NSCPDFs.totalIPDFVals = NSCPDFs.totalIPDFVals(NSCInds);
        NSCPDFs.phot1PDFVals = NSCPDFs.phot1PDFVals(NSCInds);
        NSCPDFs.phot2PDFVals = NSCPDFs.phot2PDFVals(NSCInds);    
        histNBins = length(NSCPDFs.binCents);
        binWidth = NSCPDFs.binCents(2) - NSCPDFs.binCents(1);
        minXVal = min(NSCPDFs.binCents)-binWidth/2;
        maxXVal = max(NSCPDFs.binCents)+binWidth/2;
    else
        minXVal = restrictedFittingDomain(1);
        maxXVal = restrictedFittingDomain(2);
    end  
end

 


if gpuMode
    curMCFunc = @glintMCFunc_GPU;
    parallel.gpu.rng('shuffle', 'Philox4x32-10'); % Philox4x32-10 is fastest
else
    curMCFunc = @glintMCFunc;
end

fitSettings=struct();
fitSettings.hopSigmas = hopSigmas;
fitSettings.relStepSigma = relStepSigma;
fitSettings.dataDir = dataDir;
fitSettings.dataFileName = dataFileName;
fitSettings.nLoops = nLoops;
fitSettings.nSamps = nSamps;

deltaIMu = measuredStats(5);
deltaISig = measuredStats(6);
% Background (dark current), separate for each channel (i.e. detectors are different)
% These are from normalised distributions as per dark = measuredDark/mean(peakEst) 
phot1DarkSig = measuredStats(1);
phot2DarkSig = measuredStats(2);
nullChannelDarkSig = measuredStats(3);

origFinDiffRelStep = finDiffRelStep;

if ~fitAll
    hopSigmas([3 4 7 8 9]) = 0;
end
            
% Manually specify histogram bins
histBinSize = (maxXVal-minXVal)/(histNBins);
histEdgesSpecify = minXVal:histBinSize:maxXVal;
binWidth = histEdgesSpecify(3) - histEdgesSpecify(2);
binCents = histEdgesSpecify(1:end-1)+binWidth/2;
% if useNSC && ~isequal(binCents,NSCPDFs.binCents)
if useNSC && abs(mean(NSCPDFs.binCents - binCents)) > eps
    disp('ERROR! When using NSC, the histogram sampling used here must be')
    disp('the same as the NSC PDF sampling.')
    return
end
    

if estimateHistErrors
    %%%%%%%%%%% MEASURED DATA UNCERTAINTIES %%%%%%%%%%
    % This is tricky - see various tries in binTest.m. Maybe need to use kernel
    % density estimate... We do want to use original errors, to
    % account for more photon noise far from the null!
    % FOR TESTING - Make up errors proportional to Poisson error
    % errorScaling = 0.01;
    % measuredHistValsErr = 1./sqrt(measuredHistVals) * errorScaling;

    % Approach: no. of observations in a bin is the sum of Bernoulli RVs where
    % p_i (for each) is the proportion of the area under N(x_i, sigma_i) within
    % the bin boundaries.
    newHistVals = zeros(1,histNBins);
    newHistValsSE = zeros(1,histNBins);
    binWidth = histEdgesSpecify(3) - histEdgesSpecify(2);
    binCents = histEdgesSpecify(1:end-1)+binWidth/2;
    %for k =20:histNBins
    for k =1:histNBins
        ll = histEdgesSpecify(k);
        ul = histEdgesSpecify(k+1);
        newHistVals(k) = sum( (measuredNullEst >= ll) & (measuredNullEst < ul) );

        allCurPs = normcdf(ul,measuredNullEst,measuredNullEstErr) - ...
            normcdf(ll,measuredNullEst,measuredNullEstErr);
        newHistValsSE(k) = ( sum( allCurPs.*(1-allCurPs) ) )/ sqrt(newHistVals(k));


    %     newHistVals(k) = sum(binornd(1, allCurPs));
    %     newHistValsSE(k) = std(binornd(1, allCurPs));
    end

    %Undefined points should really have inf errors (i.e. have no effect)
    newHistValsSE(isnan(newHistValsSE)) = inf;

    % Normalise as PDF
    PDFnormAmt = sum(newHistVals) * binWidth;
    newHistVals = newHistVals / PDFnormAmt;
    newHistValsErr = newHistValsSE / PDFnormAmt;

    measuredHistVals = newHistVals;
    measuredHistValsErr = newHistValsErr;

    if errorScalingFactor ~= 1
        measuredHistValsErr(binCents >= errorScalingCutoff) = ...
            measuredHistValsErr(binCents >= errorScalingCutoff)*errorScalingFactor;
    end

    %measuredHistValsErr(isinf(measuredHistValsErr)) = 1e24;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    [measuredHistVals, histBins, binInds]=histcounts(measuredNullEst,histEdgesSpecify, ...
        'Normalization', 'pdf');
    measuredHistValsErr = zeros(1, length(measuredHistVals));
    binWidth = histBins(3) - histBins(2);
    binCents = histBins(1:end-1)+binWidth/2;
end


%%%%%%%%%%%%%%%%% Do Fitting %%%%%%%%%%%%%%%%
guessParams = [deltaPhiMu, deltaPhiSig, deltaIMu, deltaISig, ...
  IrMu, IrSig, phot1DarkSig, phot2DarkSig, nullChannelDarkSig, astroNull];

% NB Not sure if it's maybe suposed to be sqrt(weights) - see
% http://stackoverflow.com/questions/22043739/weighted-curve-fitting-with-lsqcurvefit
% weightsVec = sqrt(1./measuredHistValsErr);


% % Just using lsqnonlin
% % Make a cost function
% costFunction = @(fittingParams) (glintMCFunc(nSamps, nLoops, ...
%     histEdgesSpecify, fittingParams) - measuredHistVals) .* weightsVec;
% 
% fitOptions = optimoptions(@lsqnonlin,'Display','iter-detailed','FinDiffRelStep',1e-2, ...
%     'TolX', 1e-9);%, ...
%     %'Algorithm', 'levenberg-marquardt');
% 
% tic
% [fittedParams, resnorm, residual] = lsqnonlin(costFunction, guessParams, [], [], fitOptions);
% toc





allFittedParams = zeros(numHops,length(guessParams));
allGuessParams = zeros(numHops,length(guessParams));
allFinDiffRelSteps = zeros(numHops,1);
allReducedChi2 = zeros(numHops,1);
allParamsLL = zeros(numHops, length(guessParams));
allParamsUL = zeros(numHops, length(guessParams));
allFittedDStruct = struct();
allFailedIters = zeros(numHops,1);
allFittedCurves = zeros(numHops,histNBins);

figure(figNum)
clf()
overallTic = tic;

for curHop = 1:numHops
    
    disp(['---- Hop ' num2str(curHop) ' of ' num2str(numHops) ' ----'])
    abandonThisIter = false;
    
    % Use fitChiSquare, which is a wrapper for lsqnonline which easily works
    % out paramter errors, etc.

    % NB Original guessParams = [deltaPhiMu, deltaPhiSig, deltaIMu, deltaISig, ...
    %    IrMu, IrSig, phot1DarkSig, phot2DarkSig, nullChannelDarkSig, astroNull];
    % and now deltaIMu, deltaISig, phot1/2DarkSig and nullChannelDarkSig are
    % fixed

    if fitAll
        modelFun = @glintMCFuncWrapper;
        modelGuessParams = guessParams;
        modelGuessParamsAll = modelGuessParams;
    else
        modelFun = @glintMCFuncWrapper5p;
        modelGuessParams = [deltaPhiMu, deltaPhiSig, ...
            IrMu, IrSig, astroNull];
        modelGuessParamsAll = [modelGuessParams(1) modelGuessParams(2) guessParams(3) guessParams(4) ...
            modelGuessParams(3) modelGuessParams(4) guessParams(7) guessParams(8) ...
            guessParams(9) modelGuessParams(5)];
    end

    if curHop > 1
        modelGuessParamsAll = modelGuessParamsAll + ...
            randn(1,length(modelGuessParamsAll)).*hopSigmas;
        finDiffRelStep = origFinDiffRelStep * 10^(randn()*relStepSigma);
    end


    
    
    % Plot data aand guess
    figure(figNum)
    %clf()
    hold on
    errorbar(binCents, measuredHistVals, measuredHistValsErr,'r-x')
    % Plot histogram saved from readData as a double-check:
    %errorbar(nullEstBinCents, nullEstPDFVals, nullEstPDFValsErr,'c-+')
    
    

    % Call glintMCFunc to get a simulated distribution and overplot it
%     dataParams = [deltaPhiMu, deltaPhiSig, deltaIMu, deltaISig, ...
%       IrMu, IrSig, phot1DarkSig, phot2DarkSig, nullChannelDarkSig, astroNull];
    tic
    [ MCoutputPDF, MCoutputPDFXs ] = curMCFunc( nSamps, nLoops, histEdgesSpecify, ...
        modelGuessParamsAll);
    toc
    plot(MCoutputPDFXs, MCoutputPDF,'m-x')
    %hold off
    
    
    
    
    dxData = zeros(1,length(binCents));
    % fitOptions = optimoptions(@lsqnonlin,'Display','iter-detailed');%,'FinDiffRelStep',1e-2, ...
    %    'TolX', 1e-9, ...
    %    'Algorithm', 'levenberg-marquardt');
    fitOptions = optimset('Display','iter','FinDiffRelStep',finDiffRelStep);
    fitOptions.ErrorsUnknown = ignoreErrors;
    fitOptions.FitUncertainty = fitUncertainties;
    fitOptions.MaxFunEvals = 200;
    %fitOptions.Scale = 1;


    % fitOptions.LowerBound = [guessParams(1) -inf -inf -inf -inf -inf -inf -inf -inf -inf];
    % fitOptions.UpperBound = [guessParams(1) inf inf inf inf inf inf inf inf inf];

    if performFit
        
        tic
        try
            [fittedModelParams,fittedDParams,gof,stddev] = ...
                   fitChiSquare(binCents,measuredHistVals,modelFun,modelGuessParams, ...
                   dxData,measuredHistValsErr,fitOptions);
        catch
            % A bit of a kludge, but handle this by setting the
            % allFailedIters flag (to mark results from this iter to be
            % ignored later) and then setting results to 0.
            disp('------ Fitting generated an exception, abandoning this iteration ------')
            allFailedIters(curHop) = 1;   
            fittedModelParams = zeros(size(modelGuessParams));
            st = struct('du',0,'dl',0');
            if fitAll
                fittedDParams = [st st st st st st st st st st];
            else
                fittedDParams = [st st st st st];
            end
            gof = struct('chi2',0,'dof',0);
        end
        toc

        if fitAll
            fittedParams = fittedModelParams;
            curULs = fittedDParams.du;
            curLLs = fittedDParams.dl;
        else
            fittedParams = [fittedModelParams(1) fittedModelParams(2) guessParams(3) guessParams(4) ...
                fittedModelParams(3) fittedModelParams(4) guessParams(7) guessParams(8) ...
                guessParams(9) fittedModelParams(5)];
            
            curULs = [fittedDParams(1).du fittedDParams(2).du 0 0 ...
                fittedDParams(3).du fittedDParams(4).du 0 0 ...
                0 fittedDParams(5).du];
            
            curLLs = [fittedDParams(1).dl fittedDParams(2).dl 0 0 ...
                fittedDParams(3).dl fittedDParams(4).dl 0 0 ...
                0 fittedDParams(5).dl];
            
        end

    else
        fittedParams = guessParams;
        gof.chi2=0;
        gof.dof=1;
        curLLs=0;
        curULs=0;
        fittedDParams=0;
    end
    
    
    hold on
    disp('Timing glintMCFunc evaluation...')
    profile on
    tic
    fittedCurve = curMCFunc( nSamps, nLoops, histEdgesSpecify, ...
        fittedParams);
    tt = toc;
    profile off
    disp(['It took ' num2str(tt) ' seconds.'])
    %%%plot(binCents, fittedCurve,'b-+')
    plot(binCents, fittedCurve,'b-+')
    % plot(binCents, residual+measuredHistVals,'g')
    legend('Measured', 'Guess', 'Fitted')

    %plot(binCents, newHistVals, 'ms')
    %errorbar(binCents, newHistVals, newHistValsErr,'m-s')
    
    axis([-inf inf histPlotYlim(1) histPlotYlim(2)])
    hold off

    allFittedCurves(curHop, :) = fittedCurve;
    
    % guessParams
    % modelGuessParams
    % fittedParams

    disp(['Current finDiffRelStep: ' num2str(finDiffRelStep)])
    disp(['Startvals: ' num2str(modelGuessParamsAll, '%f')])
    disp(['Endvals:   ' num2str(fittedParams, '%f')])

    allFittedParams(curHop, :) = fittedParams;
    allGuessParams(curHop, :) = modelGuessParamsAll;
    allReducedChi2(curHop) = gof.chi2/gof.dof;
    allParamsLL(curHop, :) = curLLs;
    allParamsUL(curHop, :) = curULs;
    allFinDiffRelSteps(curHop) = finDiffRelStep;
    
    curFieldName = ['hop' num2str(curHop)];
    allFittedDStruct.(curFieldName) = fittedDParams;
    
    save(fitSaveFileName, 'allFittedParams', 'allGuessParams', 'allFinDiffRelSteps', ...
        'allReducedChi2', 'allParamsLL', 'allParamsUL', 'allFittedDStruct', ...
        'binCents', 'measuredHistVals', 'measuredHistValsErr', ...
        'allFailedIters', 'fitSettings', 'histEdgesSpecify', 'binWidth', ...
        'allFittedCurves');
    
    disp(['Reduced chi2: ' num2str(allReducedChi2(curHop), '%5.3g')])
    
    avIterTime = toc(overallTic)/curHop;
    disp(['Average time per hop: ' num2str(avIterTime/60) ' minutes'])
    pause(2)
end

tt=toc(overallTic);
disp(['Total execution time: ' num2str(tt/60/60) ' hours'])





