function [ outputPDF, outputPDFXs ] = glintMCFunc_GPU( nSamps, nLoops, histEdgesSpecify, ...
    dataParams)
%glintMCFunc A function (e.g. to be called by a fitting algorithm) which
% produces a PDF of simulated GLINT data via MC method.
%
% Implemented on GPU for speed!
% Uses single precisioon for speed/
%
% dataParamas = [deltaPhiMu, deltaPhiSig, deltaIMu, deltaISig, ...
%   IrMu, IrSig, phot1DarkSig, phot2DarkSig, nullChannelDarkSig, astroNull]
% 

%disp('GPU func')
%parallel.gpu.rng('shuffle', 'Philox4x32-10'); % Philox4x32-10 is fastest
dType = 'single';

histNBins = single(length(histEdgesSpecify) - 1);
% accumPdf = zeros(nLoops,histNBins, dType); 
% allPdfXs = zeros(nLoops,histNBins, dType); 
accumPdf = zeros(nLoops,histNBins, dType, 'gpuArray'); 
allPdfXs = zeros(nLoops,histNBins, dType, 'gpuArray');

dataParams = single(dataParams);
deltaPhiMu = dataParams(1);
deltaPhiSig = dataParams(2);
deltaIMu = dataParams(3);
deltaISig = dataParams(4);
IrMu = dataParams(5);
IrSig = dataParams(6);
phot1DarkSig = dataParams(7);
phot2DarkSig = dataParams(8);
nullChannelDarkSig = dataParams(9);
astroNull = dataParams(10);

%disp(dataParams)

for k = 1:nLoops
    
    IAv = single(0.5); %Average input beam intensity, keep as 0.5 for now.
    
    %%%dISamps = randn(1,nSamps)*deltaISig + deltaIMu; 
    % I think no need to use arrayfun..?
    % Also doing as column major as supposedly faster...
    %dISamps = randn(nSamps,1, dType, 'gpuArray') * deltaISig + deltaIMu; 
    r = randn(nSamps,1, dType, 'gpuArray');
    dISamps = r * deltaISig + deltaIMu;
    
    % Don't let flux go negative...
    dISamps = dISamps((dISamps >= -1) & (dISamps <= 1));
%     goodInds = (dISamps >= -1) & (dISamps <= 1);    
%     dISamps = dISamps(goodInds);
    
    actualNSamps = length(dISamps);
    dPhiSamps = randn(actualNSamps,1, dType, 'gpuArray') * deltaPhiSig + deltaPhiMu;
    IrSamps = randn(actualNSamps,1, dType, 'gpuArray')*IrSig + IrMu;
      
    nullChannelDarkSamps = randn(actualNSamps,1, dType, 'gpuArray')*nullChannelDarkSig;
    phot1DarkSamps = randn(actualNSamps,1, dType, 'gpuArray')*phot1DarkSig;
    phot2DarkSamps = randn(actualNSamps,1, dType, 'gpuArray')*phot2DarkSig;
    
    I1 = IAv*(dISamps+1);
    I2 = IAv*(1-dISamps);
    
    % Normalise dark distributions with peakest (which should be ~1).
    % I think you can do this here since mean of dark distributions is zero
    %%%IPlusNormFactor = mean(0.5*(I1 + I2 +2*sqrt(I1.*I2)));
    rvec = arrayfun(@(I1, I2) 0.5*(I1 + I2 +2*sqrt(I1.*I2)), I1, I2);
    IPlusNormFactor = mean(rvec);

    nullChannelDarkSamps = nullChannelDarkSamps/IPlusNormFactor;
    phot1DarkSamps = phot1DarkSamps/IPlusNormFactor;
    phot2DarkSamps = phot2DarkSamps/IPlusNormFactor;
    
%     dPhiSamps = dPhiSamps(goodInds);
%     IrSamps = IrSamps(goodInds);
%     nullChannelDarkSamps = nullChannelDarkSamps(goodInds);
%     phot1DarkSamps = phot1DarkSamps(goodInds);
%     phot2DarkSamps = phot2DarkSamps(goodInds);
    
    % Use the actual noise measured on ch1:
    %%%IMinusSamps = 0.5*(I1 + I2 -2*cos(dPhiSamps).*sqrt(I1.*I2)) + nullChannelDarkSamps;
    IMinusSamps = arrayfun(@(I1, I2, dPhiSamps, nullChannelDarkSamps) ...
        0.5*(I1 + I2 -2*cos(dPhiSamps).*sqrt(I1.*I2)) + nullChannelDarkSamps, ...
        I1, I2, dPhiSamps, nullChannelDarkSamps);
    
    % Now add noise to phot chans and then derive peakEst, since this is
    % what is done in actual data reduction
    I1 = I1 + phot1DarkSamps;
    I2 = I2 + phot2DarkSamps;
    %%I1 = arrayfun(@(I1, phot1DarkSamps) I1 + phot1DarkSamps, I1, phot1DarkSamps);
    %%I2 = arrayfun(@(I2, phot2DarkSamps) I2 + phot2DarkSamps, I2, phot2DarkSamps);
    
    % This is a bit of a hack, but needed to stop IPlus going complex,
    % since I can be negative due to added noise. But in reality, the noise
    % component doesn't have this sqrt relationshoip since it is not
    % coherent...
    % TODO - again, these indexing operations probably need arrayfun?
    I1(I1 < 0) = 0;
    I2(I2 < 0) = 0;
    
    %IPlusSamps = 0.5*(I1 + I2 +2*cos(dPhiSamps).*sqrt(I1.*I2));
    %%%IPlusNoDPhiSamps = 0.5*(I1 + I2 +2*sqrt(I1.*I2));
    IPlusNoDPhiSamps = arrayfun(@(I1, I2) 0.5*(I1 + I2 +2*sqrt(I1.*I2)), ...
        I1, I2);
    
    NullSamps = IMinusSamps ./IPlusNoDPhiSamps;
    %%NullSamps = arrayfun(@(a, b) a./b, IMinusSamps, IPlusNoDPhiSamps);
    
    % Add in astronomical null
    NullSamps = NullSamps + astroNull;
    
    % Multiply *null* by Ir
    NullSamps = NullSamps .* IrSamps;
    %%NullSamps = arrayfun(@(a, b) a.*b, NullSamps, IrSamps);

    % Measure pdf
    % NB histcounts supports gpuArrays
    [histVals, histBins]=histcounts(NullSamps,histEdgesSpecify,'Normalization','pdf');
    %[histVals, histBins]=histcounts(gather(NullSamps),histEdgesSpecify,'Normalization','pdf');
    binWidth = histBins(3) - histBins(2);
    binCents = histBins(1:end-1)+binWidth/2;
    nullSampMCPdfIter = histVals;
    
    allPdfXs(k,:) = binCents;
    accumPdf(k,:) = nullSampMCPdfIter;
    
end

% outputPDF = mean(accumPdf,1);
% outputPDFXs = allPdfXs(1,:);
outputPDF = gather(mean(accumPdf,1));
outputPDFXs = gather(allPdfXs(1,:));

%disp('')









