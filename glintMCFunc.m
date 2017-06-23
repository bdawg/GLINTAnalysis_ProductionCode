function [ outputPDF, outputPDFXs ] = glintMCFunc( nSamps, nLoops, histEdgesSpecify, ...
    dataParams)
%glintMCFunc A function (e.g. to be called by a fitting algorithm) which
% produces a PDF of simulated GLINT data via MC method.
%
% dataParamas = [deltaPhiMu, deltaPhiSig, deltaIMu, deltaISig, ...
%   IrMu, IrSig, phot1DarkSig, phot2DarkSig, nullChannelDarkSig, astroNull]
% 

rng('shuffle');
histNBins = length(histEdgesSpecify) - 1;
accumPdf = zeros(nLoops,histNBins); 
allPdfXs = zeros(nLoops,histNBins); 

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

parfor k = 1:nLoops
%for k = 1:nLoops
    
    IAv = 0.5; %Average input beam intensity, keep as 0.5 for now.
    dISamps = randn(1,nSamps)*deltaISig + deltaIMu;
    dPhiSamps = randn(1,nSamps)*deltaPhiSig + deltaPhiMu;
    IrSamps = randn(1,nSamps)*IrSig + IrMu;
      
    nullChannelDarkSamps = randn(1,nSamps)*nullChannelDarkSig;
    phot1DarkSamps = randn(1,nSamps)*phot1DarkSig;
    phot2DarkSamps = randn(1,nSamps)*phot2DarkSig;
    
    % Don't let flux go negative...
    goodInds = (dISamps >= -1) & (dISamps <= 1);    
    dISamps = dISamps(goodInds);
    I1 = IAv*(dISamps+1);
    I2 = IAv*(1-dISamps);
    
    % Normalise dark distributions with peakest (which should be ~1).
    % I think you can do this here since mean of dark distributions is zero
    IPlusNormFactor = mean(0.5*(I1 + I2 +2*sqrt(I1.*I2)));
    nullChannelDarkSamps = nullChannelDarkSamps/IPlusNormFactor;
    phot1DarkSamps = phot1DarkSamps/IPlusNormFactor;
    phot2DarkSamps = phot2DarkSamps/IPlusNormFactor;
    
    dPhiSamps = dPhiSamps(goodInds);
    IrSamps = IrSamps(goodInds);
    nullChannelDarkSamps = nullChannelDarkSamps(goodInds);
    phot1DarkSamps = phot1DarkSamps(goodInds);
    phot2DarkSamps = phot2DarkSamps(goodInds);
    
    % Use the actual noise measured on ch1:
    IMinusSamps = 0.5*(I1 + I2 -2*cos(dPhiSamps).*sqrt(I1.*I2)) + nullChannelDarkSamps;
    % Now add noise to phot chans and then derive peakEst, since this is
    % what is done in actual data reduction
    I1 = I1 + phot1DarkSamps;
    I2 = I2 + phot2DarkSamps;
    
    % This is a bit of a hack, but needed to stop IPlus going complex,
    % since I can be negative due to added noise. But in reality, the noise
    % component doesn't have this sqrt relationshoip since it is not
    % coherent...
    I1(I1 < 0) = 0;
    I2(I2 < 0) = 0;
    
    %IPlusSamps = 0.5*(I1 + I2 +2*cos(dPhiSamps).*sqrt(I1.*I2));
    IPlusNoDPhiSamps = 0.5*(I1 + I2 +2*sqrt(I1.*I2));

    NullSamps = IMinusSamps ./IPlusNoDPhiSamps;

    % Add in astronomical null
    NullSamps = NullSamps + astroNull;
    
    % Multiply *null* by Ir
    NullSamps = NullSamps .* IrSamps;

    % Measure pdf
    [histVals, histBins]=histcounts(NullSamps,histEdgesSpecify,'Normalization','pdf');
    binWidth = histBins(3) - histBins(2);
    binCents = histBins(1:end-1)+binWidth/2;
    nullSampMCPdfIter = histVals;
    
    allPdfXs(k,:) = binCents;
    accumPdf(k,:) = nullSampMCPdfIter;
    
end

outputPDF = mean(accumPdf,1);
outputPDFXs = allPdfXs(1,:);








