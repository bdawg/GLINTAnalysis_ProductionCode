function [ binnedData, binnedDataErr ] = binDataErr( dataToBin, binSize )
%binDataErr Bin data and calculate standard error within bin

    nRows = length(dataToBin)/binSize;
    binMat = reshape(dataToBin,binSize,nRows);
    binnedData = mean(binMat,1);
    binnedDataErr = std(binMat, 0, 1)/sqrt(binSize);
    
end
