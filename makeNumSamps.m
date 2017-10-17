function [ outputSamples ] = makeNumSamps( pdf, binCents, nSamps, ...
    interpFactor, smoothAmt )
%makeNumSamps Make numerical random samples from a given PDF

binWidth = binCents(3) - binCents(2);
dType = 'single';

% Interpolating PDF onto a finer grid...
if interpFactor > 0
    newBinCents = min(binCents):binWidth/interpFactor:max(binCents)+binWidth;
    newPdf = interp1(binCents, pdf, newBinCents, 'linear', 'extrap');
    binCents = newBinCents;
    binWidth = binCents(3) - binCents(2);
    pdf = newPdf;
end

% Smoothing PDF with moving average filter...
if smoothAmt > 0
    pdf = smooth(pdf, smoothAmt);
end

pdf = pdf / sum(pdf);
cdf = cumsum(pdf);
[cdf, mask] = unique(cdf);
cdf_binCents = binCents(mask);
cdf_binCents = cdf_binCents + binWidth/2;

uniRand = rand(nSamps, 1, dType, 'gpuArray');
outputSamples = interp1(cdf, cdf_binCents, uniRand);

end

