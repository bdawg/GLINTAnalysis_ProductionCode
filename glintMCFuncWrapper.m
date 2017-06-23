function [ result ] = glintMCFuncWrapper( params, xData )
%glintMCFuncWrapper Wrapper for glintMCFunc to make it work with
%fitChiSquare
% NOTE - This actually ignores the xData passed, and instead uses the
% global histEdgesSpecify variable, from which xData is derived.

global nSamps nLoops histEdgesSpecify guessParams

% fixedParams = params;
% fixedParams(paramsToFix) = guessParams(paramsToFix);

[ outputPDF, outputPDFXs ] = glintMCFunc( nSamps, nLoops, histEdgesSpecify, ...
    params);

result = outputPDF';

end

