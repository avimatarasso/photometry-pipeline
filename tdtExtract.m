% Code written by Christian Pedersen, adapted by Avi Matarasso
% Michael Bruchas Lab - UW

% Extracts photometry data from TDT Data Tank and outputs time series as
% .mat file
% Assumes tankname (such as Avi_stim-210913-102052) is in path_to_data 
% REQUIRES tdt2mat.m file to be in same folder (filepath)

%% Reset MatLab workspace - clears all variables and the command window

clear all;  % clear all variables
close all;  % close all open graphs
clc   % clear command window

% add the path with your code and scripts
addpath('G:\code\Updated_package_102025')
% CUSTOMIZE

%% Only edit this section

% point to tanks (enter file path)
%path_to_data = 'G:\PhotomData\NEED TO ANALYZE\STIM ASSAYS\LCterminals\21060811_15Hz3sec_10Hz30sec\Avi_autoStim-210218-161112\dlight_CA1'; %must end in \
%savename = 'dlight_BLA';

path_to_data = 'G:\PhotomData\DBHKO\LC-CA1stim_GRABDA\5hz3s, 20hz3s\Temp Testing FOlder';
savename = 'GRABDA_CA1'; % this will be the 
subjectPattern = ['*' savename '*'];

% Change subON to 1 if you want to subtract the LLS fit 
% Must customize if you would like to subtract 405 signal from 470 signal 
details.subON  = 1; %1 for both sensors; if 0,then normalize and center around Sesnor 
details.subFit = 1; %1 to subtract the linear least squares fit
%if 0, then normalize to sensor and subtract iso alone

% stimON and fearON should be 1 if you want the TTL from fear box or you
% want recorded stim points, leave as zero otherwise
% are you using stim or fearTTL? if you do both, you might need to adjust code
details.stimON = 1;
details.fearON = 0;
details.twomice = 0;
details.dualColor = 0;
details.twoGFP = 0; filenameA_replace = 'CA1'; 
filenameB_replace = 'BLA'; % for 470B recording
 % Bruchas lab specifics
details.static = 0; %1 if run on static system in Bruchas lab
details.auto   = 0; %1 if run with automated stim protocol - make sure tank name matches one below

% do you want all figures to be shown? allfigures == 1 will show all
% check = 1 will show you raw figs, linear least squares fit, and df/f%
details.allfigures = 1;
details.check = 1;

%do not edit
saveON = 1; %save mat files
saveTXTON =1; %save txt files 0  = won't save

%define tankname for tdt extraction
tankname = 'Avi_autoStim-251001-140754'; %'Avi_stim-220922-102209';
if details.auto
    tankname = 'Avi_autoStim-210907-110207'; % 'Avi_autoStim-210907-110207';    
elseif details.dualColor
    tankname = 'Avi_dual_color-230125-083718';
elseif details.twoGFP
    tankname = 'Avi_TwoGFP_TwoSubjects-250818-111005';
end

% What period of time do you want to baseline to? This will be your F0
details.BLlength = 160;

%do you want to save the current script, 1 = yes
saveCurrent = 1;

%% save current script
if saveCurrent
    FileNameAndLocation = [mfilename('fullpath')];
    if ~strcmp(path_to_data(end),'\')
        path_to_data = [path_to_data '\'];
    end
    scriptName = strsplit(FileNameAndLocation,'\'); scriptName = [path_to_data scriptName{end}]; 
    currentTime = datestr(clock, 30); currentTime = [currentTime(1:8) '_' currentTime(10:13)];
    newbackup=sprintf(['%s' '_dt' currentTime '.m'],scriptName);
    currentfile=strcat(FileNameAndLocation, '.m');

    copyfile(currentfile,newbackup);
end

%% set up Tank name, variables to extract
tankdir = [path_to_data];

% go to your working directory / data path
workPath = [path_to_data '\' tankname];
cd(workPath)
workdir = dir(subjectPattern); %dir(pwd) for all directories

% Get a logical vector that tells which items are a directory in this folder.
dirFlags = [workdir.isdir];

% Extract only those that are directories.
workdir = workdir(dirFlags);
% Allows us to ignore temporary files in the directory 
workdir = workdir(~startsWith({workdir.name}, '.'));
workdir = workdir(~startsWith({workdir.name}, 'vars'));
txtfiles = dir(['*.txt'])';

for ii = 1:length(workdir)
%tdt files must match protocol
blockname = workdir(ii).name; % name of your file
if ~isempty(strfind(blockname,'.mat'))
    blockname = blockname(1:end-4); %get rid of the extension (.mat) at the end
end

% pick the file name to save time series (must end in .mat)
filename = [blockname '.mat'];

subject = strsplit(blockname,'_'); subject = subject{1};

%convert all of the files 
if saveTXTON
    [datOUT,Dts,timing] = tdtCONVERT(details, tankdir, tankname, blockname);
else
    [datOUT,Dts] = tdtCONVERT(details, tankdir, tankname, blockname);
end


%sampling Frequency
Fs = round(1/(max(Dts)/length(Dts))); %sample frequency Hz
details.Fs = Fs; 

% preprocessing with applying linear least square and exponential
if ~details.twomice &&  ~details.twoGFP
    if details.dualColor && ~details.static
        dat435 = datOUT.dat435;
        dat490 = datOUT.dat490;
        dat565 = datOUT.dat565;
        SensorStr='490'
        IsoStr ='435'
        [dataFiltGFP,fitIso,f2,fitcurve,dataFiltExp] = fitAndSub(dat490, dat435, Dts, details, SensorStr, IsoStr,subject);
        dfFGFP = dataFiltGFP;
        SensorStr='565';
        IsoStr ='435';
        [dataFiltRFP,fitIso,f2,fitcurve,dataFiltExp] = fitAndSub(dat565, dat435, Dts, details, SensorStr, IsoStr,subject);    
        dfFRFP = dataFiltRFP;
    elseif details.dualColor && details.static
        dat405 = datOUT.dat435.dat405A; 
        dat465 = datOUT.dat490.dat470A;
        dat560 = datOUT.dat565;
        SensorStr='465';
        IsoStr ='405';
        [dataFiltGFP,fitIso,f2,fitcurve,dataFiltExp] = fitAndSub(dat465, dat405, Dts, details, SensorStr, IsoStr,subject);
        dfFGFP = dataFiltGFP;
        SensorStr='560';
        IsoStr ='405';
        [dataFiltRFP,fitIso,f2,fitcurve,dataFiltExp] = fitAndSub(dat560, dat405, Dts, details, SensorStr, IsoStr,subject);    
        dfFRFP = dataFiltRFP;
    else
        dat405 = datOUT.dat405; 
        dat470 = datOUT.dat470;
        SensorStr='470';
        IsoStr ='405';
        [dataFilt,fitIso,f2,fitcurve,dataFiltExp] = fitAndSub(dat470, dat405, Dts,details, SensorStr, IsoStr,subject);
        dfF = dataFilt; %already dF/F %(dataFilt - median(dataFilt))./abs(median(dat470)); % this gives deltaF/F
        dfFPerc = dfF.*100; % make a percentage
        data1 = dfFPerc; %data1
    end

elseif details.twoGFP
    dat470A = datOUT.dat470A;
    dat470B = datOUT.dat470B;
    dat405A = datOUT.dat405A;
    dat405B = datOUT.dat405B;  

    SensorStr='470';
    IsoStr ='405';
    if length(dat405A)>length(Dts) 
        dat405A = dat405A(1:lengthDts);
        dat470A = dat470A(1:lengthDts);
    end
    [dataFiltA,fitIsoA,f2A,fitcurveA,dataFiltExpA] = fitAndSub(dat470A, dat405A, Dts, details, SensorStr, IsoStr,subject);
    
    SensorStr='470';
    IsoStr ='405';
    if length(dat405B)>length(Dts) 
        dat405B = dat405B(1:lengthDts);
        dat470B = dat470B(1:lengthDts);
    end
    [dataFiltB,fitIsoB,f2B,fitcurveB,dataFiltExpB] = fitAndSub(dat470B, dat405B, Dts, details, SensorStr, IsoStr,subject);
    
    dfFA = dataFiltA;
    dfFB = dataFiltB; %already dF/F %(dataFilt - median(dataFilt))./abs(median(dat470)); % this gives deltaF/F
    dfFPercA = dfFA.*100; % make a percentage
    dfFPercB = dfFB.*100; % ma

else
    SensorStr = '470D';
    IsoStr    ='405D';
    dat405D = datOUT.dat405.dat405D; 
    dat470D = datOUT.dat470.dat470D;
    if length(dat405D)>length(Dts) 
        dat405D = dat405D(1:lengthDts);
        dat470D = dat470D(1:lengthDts);
    end
    [dataFiltD,fitIsoD,f2D,fitcurveD,dataFiltExpD] = fitAndSub(dat470D, dat405D, Dts, details, SensorStr, IsoStr,subject);

    SensorStr = '470C';
    IsoStr    ='405C';
    dat405C = datOUT.dat405.dat405C; 
    dat470C = datOUT.dat470.dat470C;
    if length(dat405C)>length(Dts) 
        dat405C = dat405C(1:length(Dts));
        dat470C = dat470C(1:length(Dts));
    end
    [dataFiltC,fitIsoC,f2C,fitcurveD,dataFiltExpC] = fitAndSub(dat470C, dat405C, Dts, details, SensorStr, IsoStr,subject);
    dfFC = dataFiltC;
    dfFD = dataFiltD; %already dF/F %(dataFilt - median(dataFilt))./abs(median(dat470)); % this gives deltaF/F
    dfFPercC = dfFC.*100; % make a percentage
    dfFPercD = dfFD.*100; % make a percentage
end

%% Save file as .mat file with specified filename
if saveON && saveTXTON
    if details.fearON && ~details.stimON && ~details.twoGFP
        save(filename,'dat470','dat405','fitIso', 'dfF', 'dfFPerc','Dts','timing', 'f2', 'dataFiltExp');
    elseif ~details.fearON && details.stimON
        save(filename,'dat470','dat405','fitIso', 'dfF', 'dfFPerc','Dts','timing', 'f2', 'dataFiltExp');
    elseif details.fearON && details.stimON
       save(filename,'dat470','dat405','fitIso', 'dfF', 'dfFPerc','Dts','timingFear','timingStim', 'f2', 'dataFiltExp');
    elseif details.dualColor
       save(filename,'dfFGFP','dfFRFP','dat435','fitIso','Dts', 'f2', 'dataFiltExp');
    elseif details.twoGFP
        warning("using Avi's twoGFP convention")

        locs_und  = regexp(filename,'_'); 
        locs_dash = regexp(filename,'-'); 
        locs2Replace = locs_und(2)+1:locs_dash(2)-1;
        filenameB = filename; filenameB(locs2Replace) = filenameB_replace; 
        filenameA = filename; filenameA(locs2Replace) = filenameA_replace; 
        %txt to replace from top 
        dat470 = dat470A;
        dat405 = dat405A;
        fitIso = fitIsoA;
        dfF = dfFA;
        dfFPerc = dfFPercA;
        f2 = f2A;
        dataFiltExp = dataFiltExpA;
        save(filenameA,'dat470','dat405','fitIso', 'dfF','dfFPerc','Dts', 'f2', 'timing', 'dataFiltExp');
        dat470 = dat470B;
        dat405 = dat405B;
        fitIso = fitIsoB;
        dfF = dfFB;
        dfFPerc = dfFPercB;
        f2 = f2B;
        dataFiltExp = dataFiltExpB;
        save(filenameB,'dat470','dat405','fitIso', 'dfF','dfFPerc','Dts', 'f2', 'dataFiltExp');

    elseif details.twomice
        warning('using Avi"s twomice convention')
        locs=regexp(filename,'_');
        endLocs = regexp(filename,'-');
        
        %Fibpho2 is 465C and 405C = BoxD on left
        BoxD_filename = [filename([1:locs(2)-1 endLocs(1):end-4]) '.mat'];
        dat470 = dat470C;
        dat405 = dat405C;
        fitIso = fitIsoC;
        dfF = dfFC;
        dfFPerc = dfFPercC;
        f2 = f2C;
        dataFiltExp = dataFiltExpC;
        details.boxName = 'BoxD'; 
        details.AnalysisDatetime = currentTime;
        details.fileName = filename;
        save(BoxD_filename,'dat470','dat405', 'details','fitIso', 'dfF','dfFPerc','Dts', 'f2', 'dataFiltExp');
        BoxC_filename = [filename([1:locs(1) locs(2)+1:end-4]) '.mat'];
        dat470 = dat470D;
        dat405 = dat405D;
        fitIso = fitIsoD;
        dfF = dfFD;
        dfFPerc = dfFPercD;
        f2 = f2D;
        dataFiltExp = dataFiltExpD;
        details.boxName = 'BoxC'; 
        details.AnalysisDatetime = currentTime;
        details.fileName = filename;
        save(BoxC_filename,'dat470','dat405','details','fitIso', 'dfF','dfFPerc','Dts', 'f2', 'dataFiltExp');
    else        
        save(filename,'dat470','dat405','details','fitIso', 'dfF','dfFPerc','Dts', 'f2', 'dataFiltExp');
    end  
elseif saveON
    if details.fearON && ~details.stimON
        save(filename,'dat470','dat405','details','fitIso', 'dfF', 'dfFPerc','Dts', 'f2', 'dataFiltExp');
    elseif ~details.fearON && details.stimON && details.subON
        save(filename,'dat470','dat405','details','fitIso', 'dfF', 'dfFPerc','Dts', 'f2', 'dataFiltExp');
    elseif ~details.fearON && details.stimON && ~details.subON
        save(filename,'dat470','dat405','details', 'dfF', 'dfFPerc','Dts');
    elseif details.fearON && details.stimON
       save(filename,'dat470','dat405','details','fitIso', 'dfF', 'dfFPerc','Dts', 'f2', 'dataFiltExp');
   elseif details.dualColor && ~details.static
       save(filename,'dfFGFP','dfFRFP','dat435','details','fitIso','Dts', 'f2', 'dataFiltExp');
   elseif details.dualColor && details.static
       save(filename,'dfFGFP','dfFRFP','dat405','details','fitIso','Dts', 'f2', 'dataFiltExp');
   elseif details.twomice
        warning('using Avi"s twomice convention')
        locs=regexp(filename,'_');
        endLocs = regexp(filename,'-');
        BoxD_filename = [filename([1:locs(2)-1 endLocs(1):end-4]) '.mat'];
        dat470 = dat470C;
        dat405 = dat405C;
        fitIso = fitIsoC;
        dfF = dfFC;
        dfFPerc = dfFPercC;
        f2 = f2C;
        dataFiltExp = dataFiltExpC;
        save(BoxD_filename,'dat470','dat405','details','fitIso', 'dfF','dfFPerc','Dts', 'f2', 'dataFiltExp');
        BoxC_filename = [filename([1:locs(1) locs(2)+1:end-4]) '.mat'];
        dat470 = dat470D;
        dat405 = dat405D;
        fitIso = fitIsoD;
        dfF = dfFD;
        dfFPerc = dfFPercD;
        f2 = f2D;
        dataFiltExp = dataFiltExpD;
        save(BoxC_filename,'dat470','dat405','details','fitIso', 'dfF','dfFPerc','Dts', 'f2', 'dataFiltExp');
   else
        save(filename,'dat470','dat405','details','fitIso', 'dfF','dfFPerc','Dts', 'f2', 'dataFiltExp');
   end  
end
end
    

%% create txt files

if saveTXTON
    subStr = savename;%[subjectPattern(1:end-1) '-']; %what does your extracted subject files start with?
    matFiles = dir(['*.mat']);

    for k = 1:length(matFiles)
        clearvars -except subStr matFiles k saveON subjectPattern saveTXTON timing savename 
        
      
        
        % extract metadata from name
        tmpName   = matFiles(k).name; tmp = strsplit(tmpName,'_'); 
        %tmp2 = tmp{3}; tmp2 = strsplit(tmp2,'-'); tmpDate = tmp2{3}; tmp2 = tmp2{1};
        %subStr = tmp{1}; region = [tmp{2} '_' tmp2];
        
        locs_und  = regexp(tmpName,'_'); 
        locs_dash = regexp(tmpName,'-'); 
        locs2add = locs_und(2)+1:locs_dash(2)-1;

        %normal
        tmp2 = tmp{1}; tmp2 = strsplit(tmp2,'-'); 
        tmpDate = tmp{3}; tmpDate = strsplit(tmpDate,'-'); tmpDate = tmpDate{2}; 
        
        subStr = [tmp2{1} '-' tmp2{2}];
         region = [savename '_' tmpName(locs2add)];

        subNumb   = [tmpName(strfind(tmpName, subStr) + length(subStr))]; % will get sub
        load(tmpName)
%FOR THE ONES ELENA NAMED LOL
%subStr = tmpName(1:8);
%timingFile = [subStr '.txt'];

        timingFile = [subStr '_' region '_' tmpDate '.txt']; %[subStr 'Stim_times' subNumb '.txt']
        fid = fopen( timingFile, 'wt' );
        if exist('timingStim')
            for i = 1:length(timingStim)
                fprintf( fid, [num2str(timingStim(i)) ' ' num2str(ones(length(timingStim(i)),1))]);
                if i ~= length(timingStim)
                    fprintf(fid,'\n');
                end
            end
        elseif exist('timingFear')
            for i = 1:length(timingFear)
                fprintf( fid, [num2str(timingFear(i)) ' ' num2str(3*ones(length(timingFear(i)),1))]);
                if i ~= length(timingFear)
                    fprintf(fid,'\n');
                end
            end
        elseif exist('timing')
            for i = 1:length(timing)
                fprintf( fid, [num2str(timing(i)) ' ' num2str(3*ones(length(timing(i)),1))]);
                if i ~= length(timing)
                    fprintf(fid,'\n');
                end
            end
        end
        fclose(fid);
    end
end

%% Save All Figures
answer = questdlg('Would you like to save figures?', ...
	'Saving?', ...
	'Yes','No','No');
switch answer
    case 'Yes'
        figureSaveName = regexprep(savename, {'*'}, {''});
        if exist('varsAndFigs','dir')
            cd('varsAndFigs')
        else
            mkdir('varsAndFigs')
            cd('varsAndFigs')
        end
        tempdir = pwd;
        FolderName = tempdir;   % Your destination folder
        FigList = findobj(allchild(0), 'flat', 'Type', 'figure'); FigList2 = FigList(1:3);
        savefig(FigList2, fullfile(FolderName,[figureSaveName '.fig']));
        %save([figureSaveName '.mat']);
        disp('Figures have been saved!')        
        cd('..')
    case 'No'
        disp('You may manually save figures if you want.')
end


