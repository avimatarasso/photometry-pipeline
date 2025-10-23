% general photometry processing script for preprocessed TDT files
% preprocessed .mat files should contain "Dts", "data1" 
%  
% please use tdtExtract_clean to extract code
% 
% please email akmat@uw.edu with any questions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% add the path to all code and data directories
%codeDir = 'C:\Users\avima\OneDrive\Documents\MATLAB\Starter Code - photometry'; %CUSTOMIZE
%dataDir = 'G:\PhotomData\NEED TO ANALYZE\STIM ASSAYS\VTAterminals\20hz3s_5hz30s\Avi_autoStim-210907-110207'; %CUSTOMIZE
PhotomPack = (genpath('G:\code\Updated_package_102025'));

addpath(PhotomPack); %addpath(dataDir);

    
%% CUSTOMIZE this whole section
%do you only want to look at the averaged data of *some* of the subjects?
% if you have subjects 0,1,2,4 in the directory you're working in, 
% and you want to look at 1 and 4, type [2 4];

%paths must be in a cell array
path_to_data = {'C:\Users\bruchasadmin\Desktop\20hz30s5hz3s\Avi_autoStim-220623-085600'};
path_to_data = {'G:\PhotomData\DBHKO\LC-CA1stim_GRABDA\5hz3s, 20hz3s\Temp Testing FOlder\Avi_autoStim-251001-140754'}; %LCinhibCNO
%%%%%%%% change the savename to get each timelocked trial (usually of stim)
doYouWantAll = 1; % make 0 if you want to only choose some of the subs %and define which ones you want through whichSub
% if doing multip[le paths, want to be careful and check whichSub you choose
whichSub = [1 2]; % which subjects do you want to look at/for multiple regions include all

%%% events what event do you want to look at
eventLabel  = 3; %i usually do 2 for tail lift and 1 for stim %-2 = leaveClosedTS, 1 = loco.startTS
manualEvent = 0; %1 if the scoring was done manually, 0 if TTL

%How many trials and do you want all?
allTrials = 1; %if you don  't want all trials, make allTrials = 0
if ~allTrials  
    trialsKept   = 1; % how many of eventLabel# did you do?
    trialsOI = 1:trialsKept; % 1:10 for example, which trials are you interested in
    trialsOI = 1; % 1:10 for example, which trials are you interested in
end

%customize these lines + defineSubject.m for your specific metadata + style
%%%% commonStr == common string between your experiments of interest
region = 'GRABDA_CA1'; 
commonStr = region;
exp = '20hz'; treatment = 'HET_UPDATED_PACK';

currDate = datestr(date,'yymmdd');
savename = [region '_' exp '_' treatment '_' currDate];
a = findstr(region,'_'); b= findstr(exp,'_');
region2 = region; region2(a) = '-'; exp(b) = ' '; 
titleName = [region2 ' ' exp ' ' treatment]; %' w/o VTAinhib'];
titleName = strrep(titleName,'_',' ');
%{
region = 'GRABNE_NAc'; commonStr = region; %put in format of grabNE_CA1
%treatment = 'SAL'; 
exp = '20hz30s';
savename = [region '_' exp]
a = findstr(region,'_'); region2 = region; region2(a) = '-'; %
titleName = [region2 ' ' exp ''];
%}
if contains(region,'dual')
    FP = 'GFP';
    if strcmp(FP, 'GFP')  
        sensor = 'GRABNE';
    else
        sensor = 'GRABDA';
    end
    savename = [region '_' sensor '_' treatment '_' exp '']
    titleName = [region2 ' ' sensor ' ' exp ' ' treatment ''];
end
%regionON = 0; %mark as 0 if region not included in name
%savename = [region '_taillift' FP]; %3

%titleName = [region2 ' ' FP ' response to tail lift'];
%titleName = ['VTA stim,' region2 ' response at 20hz for 3s'];
%titleName = ['tdtmto control stim for 20hz 30s, ' region2 'response'];
%titleName = ['LC terminal stim at 20hz for 3s, dlight in vCA, SAL'];

%Do you want to adjust so its just the most recent headentry
if contains(savename,'headEntry')
    headEntryQ = 1; %make sure eventLabel is 2!
if headEntryQ
    warning('headEntryQ is 1!')
    if eventLabel ~= 2
        warning('eventLabel is not 2!')
    end
end
else
    headEntryQ = 0;
end

%%% do you want to z-score? if yes, make zQuest = 1, if df/f, make zQuest = 0
zQuest = 1; zBL = 0;
% do you want not All of the Trials? if yes ==1. If negative, will do the last ones 
normalQ = 1; %do you want to normalize to time immediately prior to event
oldQ = 0; % do you want to use the old way of analyzing data (Christian's way)

%%% do you want each subject's timelocked and averaged traces and heat map?
% if averageQuest = 1, they will be shown when you run the code
averageQuest = 0;

%%% how much time (in s) after the event do you want to look at 
%timewindow = 30; % +/- seconds from event %usually 120
%specialCase = 2; %multiplier for timewindow to end. [-timewindow specialCase*timewindow]
% If only want to see time within +/- timewindow, keep specialCase as 1
%older options above, use timeBefore and after below instead
timeBefore = 30;%2*60;% 30; %30 to 60 for short, 30 to 180
timeAfter  = 60;%7*60; 
xlims = [-timeBefore timeAfter];
delay = 0;

%CUSTOMIZE
dsFactor = 500; % how much do you want to downsample by?
baselineTime = 60; % in s
FS = 1017.25; % sampling frequency of synapse

%% Define Directory and constant variables

workdir = defineDir(path_to_data,commonStr);

%Define the working directories based on whichSub
if ~doYouWantAll
    workdir = workdir(whichSub);
end

%initialize baseline variables and timing arrays
baselineSD = zeros(length(workdir),1);
baselineMu = zeros(length(workdir),1);
subName = cell(1, length(workdir));
allTimings    = cell(length(workdir),1);
lengthOfBL = floor(abs(timeBefore)*FS);
lengthOfData = floor((timeAfter+timeBefore)*FS);
trials = []; firstTimings = []; cols = 0;

% Allows us to ignore temporary files in the directory 
workdir = workdir(~startsWith({workdir.name}, '._'));
%% save current script
FileNameAndLocation = [mfilename('fullpath')];
if ~strcmp(path_to_data(end),'\')
    path_to_data = [path_to_data '\'];
end
scriptName = strsplit(FileNameAndLocation,'\'); scriptName = strjoin([path_to_data scriptName{end}]); 
currentTime = datestr(clock, 30); currentTime = [currentTime(1:8) '_' currentTime(10:13)];
newbackup=sprintf(['%s' '_dt' currentTime '.m'],scriptName);
currentfile=strcat(FileNameAndLocation, '.m');
copyfile(currentfile,newbackup);


%% Set up events and align them

for subjectIdx = 1:length(workdir) 

    clearvars -except sess cols k itchStack workdir savename trialname FP headEntryQ sensor ...
         tempPath path_to_data sTs tT whichSub zQuest zBL baselineTime baselineSD baselineMu...
        subName subStr averageQuest firstQuest dsFactor trialsOI region notAllTrials nn minTime ...
        eventLabel allTrials timewindow specialCase normalQ lengthOfBL lengthOfData trials trialsKept oldQ ...
        timeBefore timeAfter proc time sem indData firstTimings subjectIdx manualEvent delay FS bl nData timingIdxs allData titleName xlims 
    
    %scroll through each photometry folder 
    cd(workdir(subjectIdx).folder)
    photoname = workdir(subjectIdx).name;
    load(photoname)
    
    if oldQ == 1
        data1 = dataOld;
    elseif contains(region,'dual')
        if strcmp(FP,'GFP')
            data1 = dfFGFP;
            dfF   = dfFGFP;
        else
            data1 = dfFRFP;
            dfF   = dfFRFP;
        end

        if exist('dat405','var')
            rawSensor = data1; rawIso = dat405; %static system
        else
            rawSensor = data1; rawIso = dat435; 
        end
    elseif exist('dfF','var')
        data1 = dfF; % CUSTOMIZE if you want to use a different parameter for your analysis
        rawSensor = dat470; rawIso = dat405; 
    else
        warning('dfF not calculated, probably using old data, not fit for LLS')
    end
    %CUSTOMIZE the labels and subject names you may want
    [subStr, subNumb,subjectLabel, txtName] = defineSubject(photoname, region,subjectIdx);

    if ~exist(txtName)
        txtName = dir([txtName(1:end-4) '*.txt']);
        txtName = txtName{subjectIdx}.name;
    end

    % txtName = [photoname '.txt'];
    rawtimings = dlmread(txtName); %CUSTOMIZE

    if headEntryQ
        [eventTime, numActual] = headEntryAdjust(rawtimings, eventLabel, subjectLabel, 1);
    else
         % find event times
        eventTime = rawtimings(:,1);
        eventtype = round(rawtimings(:,2));
        eventTime = eventTime(eventtype==eventLabel);
        %eventTime = eventTime(end);
    end

    if ~isempty(eventTime)
        subName{subjectIdx} = subjectLabel;
    end
    deletion_window = 20;
    eventTime=eventTime(eventTime>5); %eventTime must be after 5s for some sort of baselin
    del_idx = [0; (diff(eventTime) < deletion_window)]; %get rid of any events with less than 10 s in between
    
    % don't need to look at the timings if there is no event
    if isempty(eventTime)
        continue
    end
    eventTime = eventTime(~del_idx); % deletes events with less than 10 s in between
    
    % If you want to limit how many of the events you use change the if statement below!
    if ~allTrials
        nn = trialsKept;
        if nn> length(eventTime)
            nnTmp= length(eventTime);
            eventTime = eventTime(1:nnTmp);
        elseif nn > 0
             eventTime = eventTime(1:nn);
        elseif nn<0
            eventTime = eventTime(end+nn+1:end);%add 1 so it is the last nn      
        end
    end
    
    %add delay?
    eventTime = eventTime + delay;

    % if manually scored and using synapse, you want to adjust for
    % nonconstant frame rate
    if manualEvent
        correctTime = max(Dts);
        vidDir = photoname(1:end-4); 
        cd(vidDir); 
        vidName = dir(['Avi*' vidDir '*.avi']); vD = vidName.name;
        if strcmp(vD(1:2), '._')
            vD = vidName(2).name;
        end
        vid = VideoReader(vD); vidTime = vid.Duration;
        vidRatio = correctTime/vidTime;
        cd('..')
        allTimings{subjectIdx} = eventTime*vidRatio; %CUSTOMIZE
    else
        allTimings{subjectIdx} = eventTime;
    end

    
    
    %find baseline for each session, between baselineTime 
    baselineMu(subjectIdx) = mean(data1(ceil(eventTime(1)*FS)-(ceil(baselineTime*FS)):ceil(eventTime(1)*FS)));
    baselineSD(subjectIdx) = std(data1(ceil(eventTime(1)*FS)-(ceil(baselineTime*FS)):ceil(eventTime(1)*FS)));%std(data1(1:ceil(baselineTime*FS)));
    
    % Z score photometry data
    if zQuest
        data1 = (data1-mean(data1))./std(data1);
        warning('You z-scored your data to the full session')
    elseif zBL
        data1 = (data1-baselineMu(subjectIdx))./baselineSD(subjectIdx);
        warning('You z-scored your data to the baseline')
    else
        data1 = (data1-mean(data1))./abs(mean(data1));
        warning('df/f, NOT z-scored')
    end    
    
    %Align the data to events 
    [nData, timingIdxs, sDiff] = alignEvent(data1, FS, eventTime, timeBefore, timeAfter);

    time = linspace(1/FS,length(data1)/FS,floor(length(data1)));
    numActualTrials  = size(nData,2);
    
    if allTrials
        trials   = [trials numActualTrials]; % how many of eventLabel# did you do?
        trialsOI = 1:numActualTrials; % 1:10 for example, which trials are you interested in
    else
        numActualTrials  = trialsKept;
        trials   = [trials numActualTrials]; % how many of eventLabel# did you do?
        trialsOI = 1:nn; % 1:10 for example, which trials are you interested in
    end
    
    if numActualTrials <= length(trialsOI) 
        %cols = (subjectIdx-1)*length(trialsOI)+trialsOI;
        cols = cols(end)+trialsOI;
        cols = cols(1:numActualTrials);
    else
        cumTrial = cumsum(trials);
        if subjectIdx > 1 
            cols = cumTrial(subjectIdx-1) + trialsOI;
        else
            cols = trialsOI;
        end
    end
   

    if averageQuest
        heatFig = 200+subjectIdx;
        avgFig = 300+subjectIdx;
        plotHeatAndIndAvg(nData, heatFig, avgFig, timeBefore, timeAfter, subStr, subjectIdx, dsFactor,FS)
    end
    
    if normalQ
        [nData] = centAndNormData(nData, FS, timeBefore); %sDiff, eventTime); %normalized to pre-event
        %NcData = cData/std(data1(1:round(180*FS))); 
    end
    
    % save your data by each subject for plotting 
    allData(:,cols) = nData(:,trialsOI); 
    
    if exist('indData','var')
        minTime = floor(min(length(indData),length(data1)))/FS; 
        
        %
        if length(indData)>length(data1)
            data1PAD = padarray(data1, [length(indData)-length(data1) 0],'post');
            indData(:,subjectIdx) =  data1PAD;
        else
            indData = padarray(indData, [length(data1)-length(indData) 0],'post');            
            indData(:,subjectIdx) =  data1;
        end
        %}
    else
        indData(:,subjectIdx) =  data1;
        minTime = floor(min(length(indData),length(data1)))/FS; 
    end
    firstTimings = [firstTimings eventTime(1)];
    
    raw.rawIso(subjectIdx,:) = rawIso;
    raw.rawSensor(subjectIdx,:) = rawSensor;
    if exist('fitIso','var')
        raw.fit405(subjectIdx,:) = fitIso;
    else
        raw.fit405 = 'NOT SUBTRACTED';
    end
    raw.dfF(subjectIdx,:)    = dfF;
    raw.subName{subjectIdx} = subStr;
    
    if ~exist('sensor','var')
        sensor = region;
    end
    %plot individual averages and heat map
    figureN = subjectIdx+100;
    plotIndividual(figureN, eventTime, raw, data1, subName, subjectIdx, dsFactor, time)%(figureN, eventTime, raw, data1, nData, subName, subjectIdx, dsFactor, time,FS,xlims, region, sensor)
    
end


%% plots
%zeroidx=find(mean(allData,1)==0);
%allData(:,zeroidx) = [];

dataForPlots.allData = allData;
dataForPlots.indData = indData;
dataForPlots.minTime = minTime;
dataForPlots.firstEvent = firstTimings;
dataForPlots.timeBefore = timeBefore;
dataForPlots.timeAfter = timeAfter;
dataForPlots.dsFactor = dsFactor;
dataForPlots.FS = FS;
dataForPlots.savename = savename;
dataForPlots.titleName = titleName;
dataForPlots.workdir = workdir;
dataForPlots.trials = trials;
dataForPlots.trialsOI = trialsOI;
dataForPlots.subName= subName(~cellfun('isempty',subName));
dataForPlots.region = region;
dataForPlots.sensor = sensor;

if exist('xlims','var')
dataForPlots.xlims = xlims;
end

% create plots, might want to CUSTOMIZE
createPlots(dataForPlots)


%% save variables in folder
timevec=linspace(-timeBefore,timeAfter,length(allData));
if isfolder([region '_' sensor])
    cd([region '_' sensor])
    if exist(savename,'var')
        save(savename, 'allData','timingIdxs','timevec','path_to_data','dataForPlots','-append')
    else
        save(savename, 'allData','timingIdxs','timevec','path_to_data','dataForPlots')
    end
    cd('..')
else
    mkdir([region '_' sensor])
    cd([region '_' sensor])
    if exist(savename,'var')
        save(savename, 'allData','timingIdxs','timevec','path_to_data','dataForPlots','-append')
    else
        save(savename, 'allData','timingIdxs','timevec','path_to_data','path_to_data','dataForPlots')
    end
    cd('..')
    % close all
end
%}

%% Save All Figures
answer = questdlg('All graphs are saved. Would you like to save the individual figures too?', ...
	'Saving?', ...
	'Yes','No','No');
switch answer
    case 'Yes'
        figureSaveName = savename;
        if exist([region '_' sensor],'dir')
            cd([region '_' sensor])
        else
            mkdir([region '_' sensor])
            cd([region '_' sensor])
        end
        tempdir = pwd;
        FolderName = tempdir;   % Your destination folder
        FigList = findobj(allchild(0), 'flat', 'Type', 'figure'); FigList2 = FigList(1:end);
        savefig(FigList2, fullfile(FolderName,[figureSaveName '.fig']));
       % save([savename '.mat']);
        disp('Figures have been saved!')        
        cd('..')
    case 'No'
        disp('You may manually save figures if you want.')
end

