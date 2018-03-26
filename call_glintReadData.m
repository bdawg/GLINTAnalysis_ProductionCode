clear all
% This is just a scratch file to call glintReadData.m so you can have all
% your filenames, etc. neatly stored away here.
%
% Make sure any filenames in glintReadData.m are removed!


% filePath = '/Users/bnorris/DontBackup/GLINTdata/201608Subaru_VegaSubset/'
% % filePath = '/Volumes/pendragon1/snert/tempFromRDN/GLINT_201608_Subaru/';
% % Vega 20160812
% startTimeString1 = '20160812T221833';
% endTimeString1 = '20160812T225417'; %Everything
% startTimeString1 = '20160812T221914'; %Nice subset
% endTimeString1 = '20160812T223026'; %Nice subset

% filePath = '/Volumes/BN_DATA_OSX/NullerData/20160322_Labtests/';
% useCustomFilenameFormat = true;
% customIDString = 'acqdataLab4';
% customPrefixLength = 12;
% startTimeString1 = '20150323T001134';
% endTimeString1 = '20170323T001648';


filePath = '/Volumes/D/data_Nuller/';
filePath = 'D:\data_Nuller\';

% Good Broadband, no wiggle
startTimeString1 = '20171129T122537';
endTimeString1 = '20171129T122825';

% % Good Broadband, RandomX wiggle ampl 0.002
% startTimeString1 = '20171129T125155';
% endTimeString1 = '20171129T125443';
% 
% % Good Broadband, RandomX wiggle ampl 0.005 - get 2nd peak
% startTimeString1 = '20171129T125557';
% endTimeString1 = '20171129T130011';

% % Good Broadband, RandomX wiggle ampl 0.002 + Yampl0.001
% startTimeString1 = '20171129T130119';
% endTimeString1 = '20171129T130407';
% 
% % Good Broadband, RandomX wiggle ampl 0.005 + Yampl0.002 - get 2nd peak
% startTimeString1 = '20171129T130541';
% endTimeString1 = '20171129T131202';
% 
% % % Good Broadband, RandomX wiggle ampl 0.0005 + Yampl0.0002
% % startTimeString1 = '20171129T131317';
% % endTimeString1 = '20171129T131605';
% 
% 
% % Broadband best coupling (with drift)
% startTimeString1 = '20171201T163439';
% endTimeString1 = '20171201T164059';
% 
% % Broadband Bad coupling 'best' null (with drift)
% startTimeString1 = '20171201T170735';
% endTimeString1 = '20171201T171023'; %First 4 files
% % endTimeString1 = '20171201T171602';
% 
% % Attempt 2 Broadband Bad coupling 'best' null (with drift)
% startTimeString1 = '20171201T173112';
% endTimeString1 = '20171201T173732'; 



% MEMS Piston wiggle -  22: 0 0 0; 30: 0.1 0 0
startTimeString1 = '20171204T133926';
endTimeString1 = '20171204T134215';

% MEMS Piston wiggle -  22: 0.1 0 0; 30: 0 0 0
startTimeString1 = '20171204T134656';
endTimeString1 = '20171204T134944';

% MEMS Tip wiggle -  22: 0 0.5 0; 30: 0 0 0
startTimeString1 = '20171204T141012';
endTimeString1 = '20171204T141412';

% % MEMS Tip wiggle -  22: 0 0 0; 30: 0 0.5 0
% startTimeString1 = '20171204T150746';
% endTimeString1 = '20171204T151034';

%  
% % MEMS Tip wiggle -  22: 0 0.5 0; 30: 0 0.5 0
% startTimeString1 = '20171204T142028';
% endTimeString1 = '20171204T142316';
% 
% % MEMS Piston wiggle -  22: 0.1 0 0; 30: 0.1 0 0
% startTimeString1 = '20171204T152113';
% endTimeString1 = '20171204T152401';

% % MEMS Wiggle -  22: 0.1 0.2 0.2; 30: 0.1 0.2 0.2
% startTimeString1 = '20171204T152951';
% endTimeString1 = '20171204T153239';

% %T/T Mirror Wiggle X,Y ampl = 0.002, .001 Random PLUS
% %MEMS Wiggle -  22: 0.1 0.1 0.1; 30: 0.1 0.1 0.1
% startTimeString1 = '20171204T155748';
% endTimeString1 = '20171204T160000';


% %MEMS Wiggle - p30 0.1 only
% startTimeString1 = '20171206T173005';
% endTimeString1 = '20171206T173252';

% %MEMS Wiggle - 22: 0 0.1 0.1 ; 30: 0.1 0.1 0.1
% startTimeString1 = '20171206T173527';
% endTimeString1 = '20171206T173815';

% %MEMS Wiggle - 22: 0 0.2 0.2 ; 30: 0.1 0.2 0.2
% startTimeString1 = '20171206T174033';
% endTimeString1 = '20171206T174321';

% %MEMS Wiggle - 22: 0 0.1 0.1 ; 30: 0.1 0.1 0.1
% startTimeString1 = '20171211T171410';
% endTimeString1 = '20171211T171658';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Final set of tests (with 3 bandpasses) on 10,11&12 Dec 2017

%%% Broadband (no filter) %%%
% 01/12/2017
% Using the t/t mirror wiggler (not MEMS)
% 4:30pm - Coupling-optimised positions, 1D scan best null, STATIC
startTimeString1 = '20171201T163009';
endTimeString1 = '20171201T163257';

% 4:34pm - 1D scan best null, Wiggler Random XAmpl 0.002, YAmpl = 0.001
startTimeString1 = '20171201T163439';
endTimeString1 = '20171201T164059';

% 5:03pm - 3D scan best null, STATIC
startTimeString1 = '20171201T170307';
endTimeString1 = '20171201T170555';

% 5:31pm - 3D scan best null, Wiggler Random XAmpl 0.002, YAmpl = 0.001
startTimeString1 = '20171201T173112';
endTimeString1 = '20171201T173732';


%%% Broadband (no filter) %%%
% 06/12/2017
% Using the MEMS mirror Wiggler
% 5:35pm - 1D scan best null, With MEMS Wiggle (22:0,0.1,0.1; 30:0.1,0.1,0.1)
startTimeString1 = '20171206T173527';
endTimeString1 = '20171206T173815';

% 5:40pm - 1D scan best null, With MEMS Wiggle (22:0,0.2,0.2; 30:0.1,0.2,0.2)
startTimeString1 = '20171206T174033';
endTimeString1 = '20171206T174321';


% %%% Broadband (no filter) %%%
% 11/12/2017
% Static at best-scanned null 3:22pm
startTimeString1 = '20171211T152213';
endTimeString1 = '20171211T152501';
% With MEMS Wiggle (22:0,0.1,0.1; 30:0.1,0.1,0.1) 3:26pm
startTimeString1 = '20171211T152646';
endTimeString1 = '20171211T153059';
% ---> TO FIT? 'BinnedData_20171211T152646-20171211T153059_binsize100_peakEstBinSize100_NSC'

% % % Static at 3D scanned null 3:56pm
% % startTimeString1 = '20171211T155611';
% % endTimeString1 = '20171211T155859';
% 
% % With MEMS Wiggle (22:0,0.1,0.1; 30:0.1,0.1,0.1) 4:00pm
% % --> Histogram looks weird... drifted?
% % startTimeString1 = '20171211T155956';
% % endTimeString1 = '20171211T160244';
% 
% 
% %%% 50nm FILTER %%%
%%%%%%%%% FAINT DATA SINCE FIBER ATTENUATED! REDID LATER %%%%%%%%%%%
% Static at best-scanned null 5:06pm
startTimeString1 = '20171211T170625';
endTimeString1 = '20171211T171245';
% With MEMS Wiggle (22:0,0.1,0.1; 30:0.1,0.1,0.1) 5:14pm
startTimeString1 = '20171211T171410';
endTimeString1 = '20171211T172819';

% Static at PHOTOM ATTENUATED TO MATCH 7:54pm
startTimeString1 = '20171211T195357';
endTimeString1 = '20171211T195645';
% With MEMS Wiggle (22:0,0.1,0.1; 30:0.1,0.1,0.1) 8:02pm
startTimeString1 = '20171211T200232';
endTimeString1 = '20171211T201726';

% Back at OPT Coupling With MEMS Wiggle (22:0,0.1,0.1; 30:0.1,0.1,0.1) 8:02pm
startTimeString1 = '20171211T202138';
endTimeString1 = '20171211T203133';
% Back at OPT Coupling, STATIC 8:33 pm
startTimeString1 = '20171211T203306';
endTimeString1 = '20171211T204054';

% 20171212
% With MEMS Wiggle (22:0,0.1,0.1; 30:0.1,0.1,0.1) 12:17pm
startTimeString1 = '20171212T121734';
endTimeString1 = '20171212T122520';
% Static at best-scanned null 12:26pm
startTimeString1 = '20171212T122646';
endTimeString1 = '20171212T123306';


%%%%%%%%% NOW BRIGHT DATA 50nm, FIBER PUT BACK %%%%%%%%%%%
% % Static at best-scanned null 1:04pm (DRIFTED LOTS)
% startTimeString1 = '20171212T130449';
% endTimeString1 = '20171212T131544';
% % With MEMS Wiggle (22:0,0.1,0.1; 30:0.1,0.1,0.1) 1:19pm (SOME DRIFT)
% startTimeString1 = '20171212T131926';
% endTimeString1 = '20171212T133459';
% % With MEMS Wiggle (22:0,0.1,0.1; 30:0.1,0.1,0.1) 1:43pm (NEW NULL OPT, SOME DRIFT)
% startTimeString1 = '20171212T134349';
% endTimeString1 = '20171212T135054';


%%%%%%%% LASER at 1550nm %%%%%%%%%
% Static at best-scanned null 4:47pm (some drift)
startTimeString1 = '20171212T164659';
endTimeString1 = '20171212T165403';
% With MEMS Wiggle (22:0,0.1,0.1; 30:0.1,0.1,0.1) 4:56pm (some drift)
startTimeString1 = '20171212T165556';
endTimeString1 = '20171212T170301';
% With SMALLER MEMS Wiggle (22:0,0.05,0.05; 30:0.05,0.05,0.05) 5:05pm
startTimeString1 = '20171212T170547';
endTimeString1 = '20171212T171042';
% With MEMS Wiggle (22:0,0.1,0.1; 30:0.1,0.1,0.1) 6:09pm (another set)
startTimeString1 = '20171212T180913';
endTimeString1 = '20171212T181533';



%%%%%%%%%%%%%%%%%%% 2016 March SCExAO Measurements %%%%%%%%%%%%%%
filePath = '/Volumes/BN_DATA_OSX/NullerData/20160322_Labtests/';
useCustomFilenameFormat = true;
customIDString = 'acqdataLab14';
customPrefixLength = 13;
startTimeString1 = '20160323T001134';
endTimeString1 = '20160323T001648';


filePath = '/Volumes/BN_DATA_OSX/NullerData/201603_OnSkyData/';
% Vega, 10 files
startTimeString1 = '20160319T053011';
endTimeString1 = '20160319T053810';

% alf Lyn - second set, shoudl be ok (10 files).
startTimeString1 = '20160318T235416';
endTimeString1 = '20160319T000216';

% alf Boo, latter set, best conditions (31 files)
startTimeString1 = '20160319T032628';
endTimeString1 = '20160319T034741';
yr = [-8e-4, 150e-4];
% startTimeString1 = '20160319T015746'; % ALL alf Boo data, lost of RM, loop open, etc.

%alf Her, good files (45 files)
startTimeString1 = '20160321T041027';
%%startTimeString1 = '20160321T043000';
endTimeString1 = '20160321T044159';
% yr = [-8e-4, 200e-4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Temp entries for convenience:

% filePath = '/Users/bnorris/DontBackup/GLINTdata/declabtest_subset/';
% With MEMS Wiggle (22:0,0.1,0.1; 30:0.1,0.1,0.1) 1:43pm (NEW NULL OPT, SOME DRIFT)
startTimeString1 = '20171212T134349';
endTimeString1 = '20171212T135054';
% Back at OPT Coupling With MEMS Wiggle (22:0,0.1,0.1; 30:0.1,0.1,0.1) 8:02pm
startTimeString1 = '20171211T202138';
endTimeString1 = '20171211T203133';

filePath = 'D:\data_Nuller\';
% With MEMS Wiggle (22:0,0.1,0.1; 30:0.1,0.1,0.1) 4:56pm (some drift)
% startTimeString1 = '20171212T165556';
% endTimeString1 = '20171212T170301';

% 4:34pm - 1D scan best null, Wiggler Random XAmpl 0.002, YAmpl = 0.001
startTimeString1 = '20171201T163439';
endTimeString1 = '20171201T164059';


 
% filePath = '/Volumes/BN_DATA_OSX/NullerData/20160322_Labtests/';
% useCustomFilenameFormat = true;
% customIDString = 'acqdataLab4';
% customPrefixLength = 12;
% startTimeString1 = '20160322T223255';
% endTimeString1 = '20160322T223906';


glintReadData;

