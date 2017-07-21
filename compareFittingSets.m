% Plot the fitted params of several different fits to compare nSamps

dataDir = 'fittedParamsFiles\';
inFiles = {'fittedParams_GPUfit_sig0_2^24its', ...
    'fittedParams_fittedParams_GPUfit_sig0_2^24itsX4', ...
    'fittedParams_fittedParams_GPUfit_sig0_2^26its', ...
    'fittedParams_fittedParams_GPUfit_sig0_2^24itsX32', ...
    'fittedParams_fittedParams_OrigStartPoints_GPUfit_sigX10_2^24itsX128', ...
    };
    %'fittedParams_wip'};

xr = [-inf, inf];
xr = [0,50];


paramsToPlot = [1 2 5 6 10];
% paramsToPlot = [1 2];
% paramsToPlot = [5];

nFiles = length(inFiles);
nPars = length(paramsToPlot);
% allNLoops = zeros(nFiles, 1);
% allNSamps = zeros(nFiles, 1);

figure(1)
clf()
for m = 1:nPars
    subplot(nPars,1,m)
    hold on
    %legend;
    for k = 1:nFiles
        load([dataDir inFiles{k}])
    %     allNLoops(k) = fitSettings.nLoops;
    %     allNSamps(k) = fitSettings.nSamps;

        nSamps = fitSettings.nLoops * fitSettings.nSamps;
        legendString = [num2str(nSamps, '%.3g') ' samples'];
        plot(allFittedParams(:,paramsToPlot(m)),'-x', 'DisplayName', legendString)
    end
    if m == 1
        legend('show')
    end
    axis([xr(1) xr(2) -inf inf])
    hold off
end

