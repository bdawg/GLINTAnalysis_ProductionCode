clear all

set(0,'DefaultAxesColorOrder', hot)

% dataDir = './fittedParamFiles/';
% %inFileName = 'fittedParams_Vega_01062016_100hops.mat';
% %inFileName = 'fittedParams_Labtest200nm_100Hops_01062016';
% inFileName = 'fittedParams_alfBoo_29052016_100hops';
% %inFileName = 'fittedParams_Labtest13_600nmWF_NoVibr_MaybeSat_runningcopy';
% inFileName = 'fittedParams_alfHer_96hops';
% inFileName = 'fittedParams_Labtest5_400nmWF_NoVibr'


dataDir = '/Users/bnorris/DontBackup/GLINTdata/20170531_Subaru/fitted/'
%inFileName = 'fittedParams_wip'
inFileName = 'fittedParams_fit5_20170531T051354-20170531T052604_binsize100_peakEstBinSize100'
inFileName = 'fittedParams_fitall_20170531T051354-20170531T052604_binsize100_peakEstBinSize100'

chisqCutoff = inf;
%chisqCutoff = 10.3;

plotBest = true; % Make a nice plot of the best fit model
plotBestErrs = false; %Plot error bars in data for the best fit plot?
forceBestIter = []; %[22]; %18,22,49,83,96 Empty array to actually use best chisq iter
showParamErrs = false;

iterRange = [1 1];
% iterRange = [1 45];
iterRange = [1 100];
%maxY = 3;

paramSelect = [1 2 5 6 10];
nullParamSelect = 10;

useSavedModelVals = true;
nLoops=16;
nSamps = 2^16;
nLoops=32;
nSamps = 2^20;
histNBins = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([dataDir inFileName]);

%%%%%%%%% TEMPORARY - is saved in future versions %%%%%%%%%%%%
% minXVal = -0.2;
% maxXVal = 1.2;
% histBinSize = (maxXVal-minXVal)/(histNBins);
% histEdgesSpecify = minXVal:histBinSize:maxXVal;
% binWidth = binCents(2) - binCents(1);


figure(1)
clf()
hold on
errorbar(binCents, measuredHistVals, measuredHistValsErr,'r-x')
%axis([min(binCents)-binWidth max(binCents)+binWidth 0 maxY])
hold off

figure(2)
clf()


goodIters = [];
allGoodParams = [];
allGoodParamsUL = [];
allGoodParamsLL = [];
allGoodChisq = [];

for k = iterRange(1):iterRange(2)
    curParams = allFittedParams(k,:);
    
    if useSavedModelVals
        curModelVals = allFittedCurves(k, :);
    else
        [ curModelVals, MCoutputPDFXs ] = glintMCFunc( nSamps, nLoops, histEdgesSpecify, ...
        curParams);
        allFittedCurves(k,:) = curModelVals;
    end
    
    curSelParams = curParams(paramSelect);
    curSelParamsLL = allParamsLL(k,paramSelect);
    curSelParamsUL = allParamsUL(k,paramSelect);
    
    %if allReducedChi2(k) <= chisqCutoff
    if true %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Disable for now, this caused problems...
        goodIters = [goodIters k];
        allGoodChisq = [allGoodChisq allReducedChi2(k)];
        allGoodParams = [allGoodParams; allFittedParams(k,:)];
        allGoodParamsUL = [allGoodParamsUL; allParamsUL(k,:)];
        allGoodParamsLL = [allGoodParamsLL; allParamsLL(k,:)];
        
        figure(1)
        hold on
        p=plot(binCents, curModelVals, 'b-');
        p.Color(4)=0.1;
        hold off
        %axis([0 1.2 0 5])
        
        disp([num2str(k) '  ' num2str(allReducedChi2(k),'%4.2f') '  : ' ...
            num2str(curSelParams,'%9.5f')])
        %disp([sprintf('%4.2f ',allReducedChi2(k)) sprintf('%9.5f ',curSelParams)])
        
        xInd = 1:length(curSelParams);
        figure(2)
        hold on
        errorbar(xInd, curSelParams, curSelParamsLL, curSelParamsUL, ...
           'x', 'MarkerSize', 12);
        hold off
        
    end
    %pause(0.1)
end

figure(1)
hold on
errorbar(binCents, measuredHistVals, measuredHistValsErr,'r-x')
hold off


[bestChisq, bestIter] = min(allGoodChisq);

disp(' ')
disp(['Best chisq: ' num2str(bestChisq) ' ---> astroNull = ' ...
    num2str(allGoodParams(bestIter,nullParamSelect),'%f') ' + ' ...
    num2str(allGoodParamsUL(bestIter,nullParamSelect),'%f') ' - ' ...
    num2str(allGoodParamsLL(bestIter,nullParamSelect),'%f')])
disp(['At iter ' num2str(bestIter)])
% disp(['Max null found: ' num2str(max(allGoodParams(:,nullParamSelect))) ...
%     ', Min: ' num2str(min(allGoodParams(:,nullParamSelect))) ...
%     ', RMS: ' num2str(std(allGoodParams(:,nullParamSelect)))])


figure(3)
clf()
hold on
[~, sortedInds] = sort(allGoodChisq);
for l = 1:length(allGoodChisq)
    errorbar(l, allGoodParams(sortedInds(l),nullParamSelect), ...
        allGoodParamsLL(sortedInds(l),nullParamSelect), ...
        allGoodParamsUL(sortedInds(l),nullParamSelect), 'o')
end
hold off


if plotBest
    if ~isempty(forceBestIter)
        plotBestIter = forceBestIter;
    else
        plotBestIter = bestIter;
    end
    
    f4=figure(14);
    clf()
    hold on
    if plotBestErrs
        errorbar(binCents, measuredHistVals, measuredHistValsErr,'r-x')
    else
        plot(binCents, measuredHistVals, 'rs')
    end
    
    bestModelVals = allFittedCurves(plotBestIter, :);
    bestPar = allGoodParams(plotBestIter,:);
    bestUL = allGoodParamsUL(plotBestIter,:);
    bestLL = allGoodParamsLL(plotBestIter,:);
    
    plot(binCents, bestModelVals, 'b-');
    
    xlabel('Null depth')
    ylabel('Normalised frequency')
    legend('Measured distribution', 'Model distribution');  
    set(gca,'FontSize',16)
    
    
    % Assemble text for annotation   
    %NB: guessParams = [deltaPhiMu, deltaPhiSig, deltaIMu, deltaISig, ...
    %       IrMu, IrSig, phot1DarkSig, phot2DarkSig, nullChannelDarkSig, astroNull];
    
    if showParamErrs
        
        fHead = '\bf Directly measured';
        fl1 = ['\rm\DeltaI\mu = ' num2str(bestPar(3), '%6.4f') ' + ' ...
            num2str(bestUL(3), '%6.4f') ' - ' num2str(bestLL(3), '%6.4f')];
        fl2 = ['\DeltaI\sigma = ' num2str(bestPar(4), '%6.4f') ' + ' ...
            num2str(bestUL(4), '%6.4f') ' - ' num2str(bestLL(4), '%6.4f')];
        fl3 = ['Phot1Dark\sigma = ' num2str(bestPar(7), '%5.3g') ' + ' ...
            num2str(bestUL(7), '%5.3g') ' - ' num2str(bestLL(7), '%5.3g')];
        fl4 = ['Phot2Dark\sigma = ' num2str(bestPar(8), '%5.3g') ' + ' ...
            num2str(bestUL(8), '%5.3g') ' - ' num2str(bestLL(8), '%5.3g')];
        fl5 = ['NullDark\sigma = ' num2str(bestPar(9), '%5.3g') ' + ' ...
            num2str(bestUL(9), '%5.3g') ' - ' num2str(bestLL(9), '%5.3g')];

        vHead = '\bf Fitted params';
        vl1 = ['\rm\Delta\phi\mu = ' num2str(bestPar(1), '%6.4f') ' + ' ...
            num2str(bestUL(1), '%6.4f') ' - ' num2str(bestLL(1), '%6.4f')];
        vl2 = ['\Delta\phi\sigma = ' num2str(bestPar(2), '%6.4f') ' + ' ...
            num2str(bestUL(2), '%6.4f') ' - ' num2str(bestLL(2), '%6.4f')];
        vl3 = ['Ir\mu = ' num2str(bestPar(5), '%6.4f') ' + ' ...
            num2str(bestUL(5), '%6.4f') ' - ' num2str(bestLL(5), '%6.4f')];
        vl4 = ['Ir\sigma = ' num2str(bestPar(6), '%6.4f') ' + ' ...
            num2str(bestUL(6), '%6.4f') ' - ' num2str(bestLL(6), '%6.4f')];
        vl5 = ['Astro null = ' num2str(bestPar(10), '%6.4f') ' + ' ...
            num2str(bestUL(10), '%6.4f') ' - ' num2str(bestLL(10), '%6.4f')];

        annotStr = {fHead, fl1, fl2, fl3, fl4, fl5, '', ...
            vHead, vl1, vl2, vl3, vl4, vl5};

        anDim = [0.48 0. 0. 0.78];
        annotation('textbox', anDim, 'String', annotStr, 'FitBoxToText', 'on', ...
            'FontSize', 12, 'Margin', 10)    
        
    else
        
        fHead = '\bf Directly measured';
        fl1 = ['\rm\DeltaI\mu = ' num2str(bestPar(3), '%6.4f')];
        fl2 = ['\DeltaI\sigma = ' num2str(bestPar(4), '%6.4f')];
        fl3 = ['Phot1Dark\sigma = ' num2str(bestPar(7), '%5.3g')];
        fl4 = ['Phot2Dark\sigma = ' num2str(bestPar(8), '%5.3g')];
        fl5 = ['NullDark\sigma = ' num2str(bestPar(9), '%5.3g')];

        vHead = '\bf Fitted params';
        vl1 = ['\rm\Delta\phi\mu = ' num2str(bestPar(1), '%6.4f')];
        vl2 = ['\Delta\phi\sigma = ' num2str(bestPar(2), '%6.4f')];
        vl3 = ['Ir\mu = ' num2str(bestPar(5), '%6.4f')];
        vl4 = ['Ir\sigma = ' num2str(bestPar(6), '%6.4f')];
        vl5 = ['Astro null = ' num2str(bestPar(10), '%6.4f')];

        annotStr = {fHead, fl1, fl2, fl3, fl4, fl5, '', ...
            vHead, vl1, vl2, vl3, vl4, vl5};

        anDim = [0.59 0. 0. 0.78];
        %anDim = [0.35 0. 0. 0.88]; %For alf Her plot
        %anDim = [0.65 0. 0. 0.82]; %For alf Lyn plot
        %anDim = [0.59 0. 0. 0.82]; %Higher
        
        annotation('textbox', anDim, 'String', annotStr, 'FitBoxToText', 'on', ...
            'FontSize', 12, 'Margin', 10)  
        
    end
    
    hold off
end









