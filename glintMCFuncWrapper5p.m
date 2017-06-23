function [ result ] = glintMCFuncWrapper5p( params, xData )
%glintMCFuncWrapper Wrapper for glintMCFunc to make it work with
%fitChiSquare
% NOTE - This actually ignores the xData passed, and instead uses the
% global histEdgesSpecify variable, from which xData is derived.

% This version only fits the 5 free params

global nSamps nLoops histEdgesSpecify guessParams

modelParams = [params(1) params(2) guessParams(3) guessParams(4) ...
    params(3) params(4) guessParams(7) guessParams(8) ...
    guessParams(9) params(5)];

[ outputPDF, outputPDFXs ] = glintMCFunc( nSamps, nLoops, histEdgesSpecify, ...
    modelParams);

result = outputPDF';

end

