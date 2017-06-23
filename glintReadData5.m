%clear all

filePrefix = 'acqdata_';
fileExtn = '.bin';

%filePath = '/data/';
%filePath = '/Volumes/silo4/snert/GLINT_Data/data/201603_OnSkyData/';
%filePath = '/Volumes/BN_DATA/NullerData/201603_OnSkyData/';
%filePath = '/Volumes/silo4/snert/GLINT_Data/2016AugustSubaru/';
filePath = '/Users/bnorris/DontBackup/GLINTdata/20170531_Subaru/';

manualDarkSpecify = [];
%manualDarkSpecify = [-9.726342e-05, -3.876065e-04, -5.243951e-04, -2.450550e-04]

darkFilesSpecify = {}; % A list of timestrings that correspond to dark files


skipRead = true; % If true, just restore the following binned file:
restoreFileName = 'BinnedData_alfBoo_ALLFILES_20160319T015746-20160319T034741_binsize100';
restoreFileName = 'BinnedData_alfBooGoodSet_20160319T032628-20160319T034741_binsize100';
restoreFileName = 'BinnedData_alfHerGood_20160321T041027-20160321T044159_binsize100';
% restoreFileName = 'BinnedData_Vega10files_20160319T053011-20160319T053810_binsize100';
%restoreFileName = 'BinnedData_alfOriAAT_20151130T023358-20151130T030113_binsize100';
% restoreFileName = 'BinnedData_Vega_WithErrs_20160319T053011-20160319T053810_binsize100_peakEstBinSize100';
restoreFileName = 'BinnedData_20170608T020709-20170608T031156_binsize100_peakEstBinSize100'


gain = 1; % V/W. Set to 1 to keep units as Volts.
subtBias = true;
nSampScl = 1; % Reduce n by this factor to account for oversampling
              % (Very rough approach - do stats properly).
              % I *think* that since std is measured after binning, we can
              % leave this at 1...
               
              
simpleNormalise = true; %For Histogram, just normalise such that max(ch) = 1
useNullEst = false; %For Histogram, use an estimated peak to estimate null depth

binSize = 100%100; % Bin together this number of samples (0 for no binning)       

binPeakEst = true; % Further bin estimates of peak (IPlus)
peakEstBinSize = 100;

doDarkHists = true; % Also do the analysis for the dark data


xr = [-inf, inf];
%xr = [0, 1.0];
%yr = [-1e-14, 10e-13]; 
yr = [-10e-4, 20e-4];
yr = [-20e-4, 100e-4];

byr = yr;
%byr = [-1e-3,1e-3]; %yrange for Binned data plots
%byr = [-inf, inf];

histNBins=100; %Number of bins for the histogram plot  
minXVal = -1.2;
minXVal = -0.2;
maxXVal = 1.2;
%histAxes = [0, 0.2];
histAxes = [-inf, inf];

doHistErrors = false; %If false, use histcounts to do it quickly.
plotDark = true; % Set to plot the signal in the dark files

saveBinnedData = true;

useCustomFilenameFormat = false;




%%% Data set:
startTimeString1 = '20160318T235416';
endTimeString1 = '20160319T000216';

%%% AAT Nov2015 - gam Ret
startTimeString1 = '20151129T234519';
%endTimeString1 = '20151129T234642'; %Just first 3 files
endTimeString1 = '20151130T004638';

%%% AAT Nov2015 - alf Ori
startTimeString1 = '20151130T023358';
endTimeString1 = '20151130T030113';


% % alf Lyn - first set, while RMs etc.
% startTimeString1 = '20160318T231025';
% endTimeString1 = '20160318T233228';
% %%%endTimeString1 = '20160318T231400'; 
% 
% % alf Lyn - second set, shoudl be ok (10 files).
startTimeString1 = '20160318T235416';
endTimeString1 = '20160319T000216';

% alf Boo, latter set, best conditions (31 files)
startTimeString1 = '20160319T032628';
endTimeString1 = '20160319T034741';
%yr = [-8e-4, 150e-4];
% %startTimeString1 = '20160319T015746'; % ALL alf Boo data, lost of RM, loop open, etc.

% % Vega, 10 files
% startTimeString1 = '20160319T053011';
% endTimeString1 = '20160319T053810';

% % eps Leo - all files, log says this is all no good
% %startTimeString1 = '20160318T203134';
% startTimeString1 = '20160318T205307'; %Starts after lab null position used
% endTimeString1 = '20160318T211216';

% %alf Her, good files (45 files)
% startTimeString1 = '20160321T041027';
% %%startTimeString1 = '20160321T043000';
% endTimeString1 = '20160321T044159';
% % yr = [-8e-4, 200e-4];

% % Altair, end of Defrere night, bad seeing.
% startTimeString1 = '20160320T055046';
% endTimeString1 = '20160320T061028';



% useCustomFilenameFormat = true;
% customIDString = 'labMeas09';
% customPrefixLength = 18;
% startTimeString1 = '20150320T221105';
% endTimeString1 = '20170320T221453';
% yr = [-1e-3, 10];
% % startTimeString1 = '20160320T222441';
% % endTimeString1 = '20160320T222611';

% %filePath = '/Volumes/BN_DATA/NullerData/20160322_Labtests/';
% filePath = '/Volumes/silo4/snert/GLINT_Data/data/20160322_Labtests/';
% useCustomFilenameFormat = true;
% customIDString = 'acqdataLab5';
% customPrefixLength = 12;
% startTimeString1 = '20150323T001134';
% endTimeString1 = '20170323T001648';
% yr = [0., 8.];
% byr = yr;

% useCustomFilenameFormat = true;
% customIDString = '2000';
% customPrefixLength = 22;
% startTimeString1 = '20160316T165951';
% endTimeString1 = '20160316T170522';
% yr = [0., 0.1];



%%%%%%%%%%%%%%%% August 2016 %%%%%%%%%%%%%%%%%

% Vega 20160814
startTimeString1 = '20160814T224045';
%endTimeString1 = '20160814T230100'; %First few files
endTimeString1 = '20160814T234820';

excludeString = [];
%excludeString = '20160323T002417'



%%%%%% beta Peg 201705
%%%yr = [-0.2e-3, 2e-3];
% Final files (after last scan), LOWFS closed:
startTimeString1 = '20170531T051354';
endTimeString1 = '20170531T052933';

% Fewer from final set:
endTimeString1 = '20170531T052604';

% Files after 2nd-last scan (still supposedly LOWFS closed)
% startTimeString1 = '20170531T045908';
% endTimeString1 = '20170531T051354';


% % Files before LOWFS closed (but at null):
% startTimeString1 = '20170531T044617';
% endTimeString1 = '20170531T045111';

excludeString = [];




%%%%%%%%%%%%%%%% AAT July 2017 %%%%%%%%%%%%%%%%%
filePath = '/Volumes/pendragon1/snert/GLINT/201706_AAT/data_Nuller/data201707_0207-0312/'
startTimeString1 = '20170608T020709';
endTimeString1 = '20170608T031156';
%endTimeString1 = '20170608T021447'; %Just first few files
endTimeString1 = '20170608T020850'; %Just first 3 files

darkFilesSpecify = {'20170608T020709', '20170608T021407', '20170608T022222' ...
    '20170608T023203', '20170608T023903', '20170608T024628', ...
    '20170608T025305', '20170608T025918', '20170608T030524', ...
    '20170608T031156'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveDir = [filePath 'BinnedSaved/'];

formatIn = 'yyyymmddTHHMMSS';
startTimeNum = datenum(startTimeString1,formatIn);
endTimeNum = datenum(endTimeString1,formatIn);

listing = dir([filePath '*' fileExtn]);
allFileNames = { listing.name };
nAllFiles = length(allFileNames);
fCount = 1;
if useCustomFilenameFormat
    pLength = customPrefixLength;
else
    pLength = length(filePrefix);
end



% Eventually put all data into 1 matrix for simplicity
nChans=4;
allData = zeros(nChans, 0);
allDataErrs = zeros(nChans, 0);
allTimes = [];
allDarkData = zeros(nChans, 0);
allDarkTimes = [];

tic
if ~skipRead
    clear ch0Av ch1Av ch2Av ch3Av ch0Averr ch1Averr ch2Averr ch3Averr
    for k = 1:nAllFiles
        skipThisFile = false;
        curName = allFileNames{k};
        %disp(curName)
        
        % If filename is in the wrong format, ignore this file
        try
            % Check if it's a dark
            dchk = strfind(curName, 'DARK');
            if length(dchk) > 0 %#ok<ISMT>
                darkFileStatus = true;
                curTimeStr = curName(pLength+6:pLength+6+15);
            else
                darkFileStatus = false;
                curTimeStr = curName(pLength+1:pLength+1+15);
            end

            % Include manually specified dark files
            if ~isempty(darkFilesSpecify)
                for s = 1:length(darkFilesSpecify)
                    check = strfind(curName, darkFilesSpecify{s});
                    if length(check) > 0 %#ok<ISMT>
                        darkFileStatus = true;
                    end
                end
            end
            
            % Check if it's the specified custom filename format
            if useCustomFilenameFormat
                cchk = strfind(curName, customIDString);
                if isempty(cchk)
                    skipThisFile = true; %Ignore file unless it matches
                end
            end
            
            % Check it doesn't match the exclude string
            if ~isempty(excludeString)
                echk = strfind(curName, excludeString);
                if ~isempty(echk)
                    skipThisFile = true;
                end
            end
            
            curTimeNum = datenum(curTimeStr,formatIn);
        catch
            disp(['Ignored file ' curName])
            skipThisFile = true;
        end
        
        
        if skipThisFile
            curTimeNum = 0; %Force this file to be skipped
        end
        
        if (curTimeNum >= startTimeNum) && (curTimeNum <= endTimeNum)
            
            isDark(fCount) = darkFileStatus
            
            % Read this file and save its statistics
            filestring=[filePath curName];
            disp(['Reading ' filestring])

            fid = fopen(filestring,'r');
            [data,count] = fread(fid,[5,inf],'double');
            fclose(fid);

            
            if length(data) > binSize*2 % Skip 'empty' files
                
                time = data(1,:);
                ch0 = data(2,:);
                ch1 = data(3,:);
                ch2 = data(4,:);
                ch3 = data(5,:);     


                % Now bin raw data and store in main array
                if binSize > 0
                    % Remove elements so bin size is evenly divisble
                    npts = length(ch0);
                    rem = mod(npts, binSize);
                    inds = 1:(npts-rem);

%                     time = downsamplebin(time(inds),2,binSize);
%                     ch0 = downsamplebin(ch0(inds),2,binSize);
%                     ch1 = downsamplebin(ch1(inds),2,binSize);
%                     ch2 = downsamplebin(ch2(inds),2,binSize);
%                     ch3 = downsamplebin(ch3(inds),2,binSize);
                    [time, dummy] = binDataErr(time(inds), binSize);
                    [ch0, ch0err] = binDataErr(ch0(inds), binSize);
                    [ch1, ch1err] = binDataErr(ch1(inds), binSize);
                    [ch2, ch2err] = binDataErr(ch2(inds), binSize);
                    [ch3, ch3err] = binDataErr(ch3(inds), binSize);

                end
                
                % Store fle averages for later plotting
                ch0Av(fCount) = mean(ch0)/gain;
                ch1Av(fCount) = mean(ch1)/gain;
                ch2Av(fCount) = mean(ch2)/gain;
                ch3Av(fCount) = mean(ch3)/gain;
                ch0Averr(fCount) = std(ch0)/sqrt(length(ch0)/nSampScl)/gain;
                ch1Averr(fCount) = std(ch1)/sqrt(length(ch1)/nSampScl)/gain;
                ch2Averr(fCount) = std(ch2)/sqrt(length(ch2)/nSampScl)/gain;
                ch3Averr(fCount) = std(ch3)/sqrt(length(ch3)/nSampScl)/gain;
                timeNums(fCount) = curTimeNum;
                

                if ~darkFileStatus
                    allData = [allData [ch0; ch1; ch2; ch3]];
                    allDataErrs = [allDataErrs [ch0err; ch1err; ch2err; ch3err]];

                    % Start new time at end of old time (so discontinuities
                    % will happen - be careful if FFTing, etc.)
                    tStep = time(2) - time(1);
                    if ~isempty(allTimes)
                        startTime = allTimes(end) + tStep;
                        time = time + startTime;
                    end
                    allTimes = [allTimes time];
                end

                if darkFileStatus
                    allDarkData = [allDarkData [ch0; ch1; ch2; ch3]];

                    % Start new time at end of old time (so discontinuities
                    % will happen - be careful if FFTing, etc.)
                    tStep = time(2) - time(1);
                    if ~isempty(allDarkTimes)
                        startTime = allDarkTimes(end) + tStep;
                        newTimes = time + startTime;
                    else
                        newTimes = time;
                    end
                    allDarkTimes = [allDarkTimes newTimes];
                end
                
                if darkFileStatus && plotDark
                    %disp(['Dark sigmas: ' num2str(std(ch0)) ', ' num2str(std(ch1)) ...
                    %    ', ' num2str(std(ch2)) ', ' num2str(std(ch3))])
                    % Plot dark data
                    savedDarkData = [ch0; ch1; ch2; ch3];
                    figure(5)
                    clf()
                    subplot(2,2,1)
                    plot(time,ch0)
                    hold on
                    plot([0 max(time)], [0 0], 'k--')
                    hold off
                    axis([xr(1) xr(2) byr(1) byr(2)]);
                    title('Channel 0 - DARK')
                    ylabel('Voltage')
                    xlabel('Time')
                    subplot(2,2,2)
                    plot(time,ch1)
                    hold on
                    plot([0 max(time)], [0 0], 'k--')
                    hold off
                    axis([xr(1) xr(2) byr(1) byr(2)]);
                    title('Channel 1 - DARK')
                    ylabel('Voltage')
                    xlabel('Time')
                    subplot(2,2,4)
                    plot(time,ch2)
                    hold on
                    plot([0 max(time)], [0 0], 'k--')
                    hold off
                    axis([xr(1) xr(2) byr(1) byr(2)]);
                    title('Channel 2 - DARK')
                    ylabel('Voltage')
                    xlabel('Time')
                    subplot(2,2,3)
                    plot(time,ch3)
                    hold on
                    plot([0 max(time)], [0 0], 'k--')
                    hold off
                    axis([xr(1) xr(2) byr(1) byr(2)]);
                    title('Channel 3 - DARK')
                    ylabel('Voltage')
                    xlabel('Time') 
                    drawnow
                end


                clear ch0;
                clear ch1;
                clear ch2;
                clear ch3;
                clear time;
                clear data;

                fCount = fCount + 1;
            else
                disp('Skipping this file - not enough data.')
            end
        end
    end  
    disp(['Total num files: ' num2str(fCount-1)])
    
    s1ix = find(isDark == false);
    s2ix = find(isDark == true);
    % Assumes set 2 is dark
    bias0 = mean(ch0Av(s2ix));
    bias1 = mean(ch1Av(s2ix));
    bias2 = mean(ch2Av(s2ix));
    bias3 = mean(ch3Av(s2ix));
    %bias = mean([bias0 bias1 bias2 bias3])
     
    if ~isempty(manualDarkSpecify) > 0
        disp('Using manually specified dark values!')
        bias0 = manualDarkSpecify(1);
        bias1 = manualDarkSpecify(2);
        bias2 = manualDarkSpecify(3);
        bias3 = manualDarkSpecify(4);
    end
    
    if subtBias
        ch0Av(s1ix) = ch0Av(s1ix) - bias0;
        ch1Av(s1ix) = ch1Av(s1ix) - bias1;
        ch2Av(s1ix) = ch2Av(s1ix) - bias2;
        ch3Av(s1ix) = ch3Av(s1ix) - bias3;
        
        % Subtract bias from dark data?
        ch0Av(s2ix) = ch0Av(s2ix) - bias0;
        ch1Av(s2ix) = ch1Av(s2ix) - bias1;
        ch2Av(s2ix) = ch2Av(s2ix) - bias2;
        ch3Av(s2ix) = ch3Av(s2ix) - bias3;
        
        allData(1,:) = allData(1,:) - bias0;
        allData(2,:) = allData(2,:) - bias1;
        allData(3,:) = allData(3,:) - bias2;
        allData(4,:) = allData(4,:) - bias3; 
        
        % Subtract bias from dark data?
        allDarkData(1,:) = allDarkData(1,:) - bias0;
        allDarkData(2,:) = allDarkData(2,:) - bias1;
        allDarkData(3,:) = allDarkData(3,:) - bias2;
        allDarkData(4,:) = allDarkData(4,:) - bias3;  
    end
end
toc

%s1ix = find(isDark == false);
%s2ix = find(isDark == true);

if skipRead
    load([saveDir restoreFileName])
end


% Now use the known chip throughout coeffs to estimate I+
% (This should be the same as ch2, but maybe less noisy ?)
ch1ToNullRatio = 1.2561;
ch4ToNullRatio = 0.9148;
phot1 = allData(1,:)/ch1ToNullRatio;
phot1(phot1 <= 0) = 1e-14; %This is actually problematic, since can be negative from noise
phot2 = allData(4,:)/ch4ToNullRatio;
phot2(phot2 <= 0) = 1e-14; %This is actually problematic, since can be negative from noise
peakEst = phot1 + phot2 +2*sqrt(phot1.*phot2); %This assumes deltaPhi=0!!


% Try binning peakEst more (it is just an estimate...)
if binPeakEst
    npts = length(peakEst);
    rem = mod(npts, peakEstBinSize);
    inds = 1:(npts-rem);
%     peakEst = downsamplebin(peakEst(inds),2,peakEstBinSize);
%     IplusTimes = downsamplebin(allTimes(inds),2,peakEstBinSize);
    
%     Note: the errors returned here ignore the errors already attached to
%     phot1 and phot2, but are re-estimated here.
    nPeakEstBins = length(peakEst(inds)) / peakEstBinSize;
    [peakEst, peakEstErr] = binDataErr(peakEst(inds),peakEstBinSize);
    [IplusTimes, ~] = binDataErr(allTimes(inds),peakEstBinSize);
    
    peakEst = interp1(IplusTimes, peakEst, allTimes,'linear','extrap');
    peakEstErr = interp1(IplusTimes, peakEstErr, allTimes,'linear','extrap');
    
    % We're adding more values to peakEstErr (via interp) but they are not 
    % really independent points, but they are correlated. Ignore this
    % correction for now (relative weights will still be right)
    %peakEstErr = peakEstErr * sqrt(binSize); 
else
    peakEstErr = zeros(1,length(peakEst));
end


%Calculate null depth estimate
nullEst = allData(2,:)./peakEst;
nullEstErr = abs(nullEst) .* sqrt( (allDataErrs(2,:)./allData(2,:)).^2 + ...
   (peakEstErr./peakEst).^2 );

% % For testing:
% % OR ignore the whole thing:
% nullEst = allData(2,:);
% simpleNormalise = true;


% Plot file statistics
figure(2)
clf()
subplot(2,2,1)
hold on
%errorbar(timeNums,ch0Av,ch0Averr)
errorbar(ch0Av(s1ix),ch0Averr(s1ix))
errorbar(ch0Av(s2ix),ch0Averr(s2ix),'r')
title('Channel 0')
xlabel('File index (arb)')
ylabel('Watts')
axis([-inf inf yr(1) yr(2)]);
hold off

subplot(2,2,2)
hold on
%errorbar(timeNums,ch1Av,ch1Averr)
errorbar(ch1Av(s1ix),ch1Averr(s1ix))
errorbar(ch1Av(s2ix),ch1Averr(s2ix),'r')
title('Channel 1')
xlabel('File index (arb)')
ylabel('Watts')
axis([-inf inf yr(1) yr(2)]);
hold off

subplot(2,2,4)
hold on
%errorbar(timeNums,ch2Av,ch2Averr)
errorbar(ch2Av(s1ix),ch2Averr(s1ix))
errorbar(ch2Av(s2ix),ch2Averr(s2ix),'r')
title('Channel 2')
xlabel('File index (arb)')
ylabel('Watts')
axis([-inf inf yr(1) yr(2)]);
hold off

subplot(2,2,3)
hold on
%errorbar(timeNums,ch3Av,ch3Averr)
errorbar(ch3Av(s1ix),ch3Averr(s1ix))
errorbar(ch3Av(s2ix),ch3Averr(s2ix),'r')
title('Channel 3')
xlabel('File index (arb)')
ylabel('Watts')
axis([-inf inf yr(1) yr(2)]);
hold off


% Print statistics
ch0AvAll = mean(ch0Av(s1ix));
ch1AvAll = mean(ch1Av(s1ix));
ch2AvAll = mean(ch2Av(s1ix));
ch3AvAll = mean(ch3Av(s1ix));
ch0AverrAll = sqrt(sum(ch0Averr(s1ix).^2)) / length(ch0Averr(s1ix));
ch1AverrAll = sqrt(sum(ch1Averr(s1ix).^2)) / length(ch1Averr(s1ix));
ch2AverrAll = sqrt(sum(ch2Averr(s1ix).^2)) / length(ch2Averr(s1ix));
ch3AverrAll = sqrt(sum(ch3Averr(s1ix).^2)) / length(ch3Averr(s1ix));
% disp(['Ch 0: ' num2str(ch0AvAll,'%e') ' +- ' num2str(ch0AverrAll,'%e')...
%     ', Bias: ' num2str(bias0, '%e')])
% disp(['Ch 1: ' num2str(ch1AvAll,'%e') ' +- ' num2str(ch1AverrAll,'%e')...
%     ', Bias: ' num2str(bias1, '%e')])
% disp(['Ch 2: ' num2str(ch2AvAll,'%e') ' +- ' num2str(ch2AverrAll,'%e')...
%     ', Bias: ' num2str(bias2, '%e')])
% disp(['Ch 3: ' num2str(ch3AvAll,'%e') ' +- ' num2str(ch3AverrAll,'%e')...
%     ', Bias: ' num2str(bias3, '%e')])
disp(['Ch 0: ' num2str(ch0AvAll,'%e') ' +- ' num2str(ch0AverrAll,'%e')])
disp(['Ch 1: ' num2str(ch1AvAll,'%e') ' +- ' num2str(ch1AverrAll,'%e')])
disp(['Ch 2: ' num2str(ch2AvAll,'%e') ' +- ' num2str(ch2AverrAll,'%e')])
disp(['Ch 3: ' num2str(ch3AvAll,'%e') ' +- ' num2str(ch3AverrAll,'%e')])


% Plot binned data
figure(3)
clf()
subplot(2,2,1)
plot(allTimes,allData(1,:))
hold on
plot([0 max(allTimes)], [0 0], 'k--')
hold off
%axis tight
axis([xr(1) xr(2) byr(1) byr(2)]);
title('Channel 0')
ylabel('Voltage')
xlabel('Time')

subplot(2,2,2)
plot(allTimes,allData(2,:))
hold on
plot([0 max(allTimes)], [0 0], 'k--')
hold off
%axis tight
axis([xr(1) xr(2) byr(1) byr(2)]);
title('Channel 1')
ylabel('Voltage')
xlabel('Time')

subplot(2,2,4)
plot(allTimes,allData(3,:))
hold on
plot([0 max(allTimes)], [0 0], 'k--')
hold off
%axis tight
axis([xr(1) xr(2) byr(1) byr(2)]);
title('Channel 2')
ylabel('Voltage')
xlabel('Time')

subplot(2,2,3)
plot(allTimes,allData(4,:))
hold on
plot([0 max(allTimes)], [0 0], 'k--')
hold off
%axis tight
axis([xr(1) xr(2) byr(1) byr(2)]);
title('Channel 3')
ylabel('Voltage')
xlabel('Time')



%%%% Now plot a histogram

% Manually specify histogram bins
histBinSize = (maxXVal-minXVal)/(histNBins);
histEdgesSpecify = minXVal:histBinSize:maxXVal;

if useNullEst
    dataForHist=nullEst;
else
    dataForHist=allData(2,:);
end

if ~doHistErrors
    figure(4)
    clf()
    if simpleNormalise
        dataForHist = dataForHist / max(dataForHist);
    end
    [histVals, histBins, binInds]=histcounts(dataForHist,histEdgesSpecify, ...
        'Normalization', 'pdf');
    binWidth = histBins(3) - histBins(2);
    binCents = histBins(1:end-1)+binWidth/2;
    axis([histAxes(1) histAxes(2) -inf inf]);
    [dummy, maxInd] = max(histVals);
    disp(['Mode: ' num2str(histBins(maxInd))])
    plot(binCents, histVals)

%     % Overplot dark histogram
%     offset = 1.05; %Arb x offset for nice plotting
%     dataForHist=savedDarkData(2,:)/mean(peakEst) + offset; %Get in about the right units
%     %dataForHist=savedDarkData(2,:);
%     hold on
%     [histVals, histBins]=histcounts(dataForHist,histEdgesSpecify);
%     binWidth = histBins(3) - histBins(2);
%     binCents = histBins(1:end-1)+binWidth/2;
%     plot(binCents, histVals*(fCount-1),'m'); %Scaled to make frequencies similar
%     [dummy, maxInd] = max(histVals);
%     hold off
    
else
    % Insert histogram with error calc here!   
end



if doDarkHists && ~skipRead
    % THIS STUFF COPIED FROM ABOVE
    phot1Dark = allDarkData(1,:)/ch1ToNullRatio;
    phot2Dark = allDarkData(4,:)/ch4ToNullRatio;
    %phot1Dark(phot1Dark <= 0) = 0; %This is actually problematic, since can be negative from noise
    %phot2Dark(phot2Dark <= 0) = 0; %This is actually problematic, since can be negative from noise
    peakEstDark = phot1Dark + phot2Dark; %No sqrt since incoherent
    %peakEstDark = phot1Dark + phot2Dark+2*sqrt(phot1Dark.*phot2Dark);
    
    if binPeakEst
        npts = length(peakEstDark);
        rem = mod(npts, peakEstBinSize);
        inds = 1:(npts-rem);
        peakEstDark = downsamplebin(peakEstDark(inds),2,peakEstBinSize);
        IplusTimesDark = downsamplebin(allDarkTimes(inds),2,peakEstBinSize);
        peakEstDark = interp1(IplusTimesDark, peakEstDark, allDarkTimes,'linear','extrap');
    end
    
    %nullEstDark=allDarkData(2,:)./peakEstDark; 
    dataForHist=peakEstDark;
    
    figure(9)
    clf()
    hold on
    if simpleNormalise
        dataForHist = dataForHist / max(dataForHist);
    end
    [histVals, histBins, binInds]=histcounts(dataForHist,histNBins, ...
        'Normalization', 'pdf');
    binWidth = histBins(3) - histBins(2);
    binCents = histBins(1:end-1)+binWidth/2;
    axis([histAxes(1) histAxes(2) -inf inf]);
    [dummy, maxInd] = max(histVals);
    %disp(['Dark Mode: ' num2str(histBins(maxInd))])
    plot(binCents, histVals)
    hold off
    
    %Normalise dark counts
    dkScaleFactor = mean(peakEst);
    
    disp(' ')
    disp(['sd of peakEstDark: ' num2str(std(dataForHist))]);
%     disp(['sd of ch0Dark: ' num2str(std(allDarkData(1,:)))]);
%     disp(['sd of ch1Dark: ' num2str(std(allDarkData(2,:)))]);
%     disp(['sd of ch2Dark: ' num2str(std(allDarkData(3,:)))]);
%     disp(['sd of ch3Dark: ' num2str(std(allDarkData(4,:)))]);
    disp(['sd of phot1: ' num2str(std(phot1Dark)/dkScaleFactor)]);
    disp(['sd of phot2: ' num2str(std(phot2Dark)/dkScaleFactor)]);
    disp(['sd of null channel: ' num2str(std(allDarkData(2,:))/dkScaleFactor)]);
    disp(['sd of bright channel: ' num2str(std(allDarkData(3,:))/dkScaleFactor)]);
    
end

% Print other required statistics
deltaIdist = (phot1 - phot2) ./ (phot1 + phot2);
disp(' ')
disp(['DeltaI mu: ' num2str(mean(deltaIdist))])
disp(['DeltaI sig: ' num2str(std(deltaIdist))])


disp(' ')
disp(['Binned sample rate (Hz): ' num2str(1/(allTimes(2)-allTimes(1)))])

% Save binned data
if saveBinnedData && ~skipRead
    measuredStats = [std(phot1Dark)/dkScaleFactor, std(phot2Dark)/dkScaleFactor, ...
        std(allDarkData(2,:))/dkScaleFactor, std(allDarkData(3,:))/dkScaleFactor, ...
        mean(deltaIdist), std(deltaIdist)];
    if binPeakEst
        peakEstFileString = ['_peakEstBinSize' num2str(peakEstBinSize)];
    else
        peakEstFileString = ''
    end
    filename = [saveDir 'BinnedData_' startTimeString1 '-' endTimeString1 '_binsize' ...
        num2str(binSize) peakEstFileString '.mat'];
    disp(['Saving binned data to ' filename])
    save(filename, 'allTimes', 'allData', 'allDataErrs', 'bias0', 'bias1', 'bias2', 'bias3', ...
        'ch0Av', 'ch1Av', 'ch2Av', 'ch3Av', 'ch0Averr', 'ch1Averr', 'ch2Averr', ...
        'ch3Averr', 's1ix', 's2ix', 'peakEst', 'peakEstErr', 'nullEst', 'nullEstErr', 'measuredStats', ...
        'histEdgesSpecify');
end


