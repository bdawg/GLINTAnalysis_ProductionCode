function [outv] = downsamplebin(invec,dimnum,amount)
    sizebefore = size(invec);
    sizemiddle = [sizebefore(1:dimnum-1),amount,sizebefore(dimnum)/amount,sizebefore(dimnum+1:length(sizebefore))];
    sizeafter = sizebefore;
    sizeafter(dimnum) = sizeafter(dimnum)/amount;
    outv = reshape(mean(reshape(invec,sizemiddle),dimnum),sizeafter);