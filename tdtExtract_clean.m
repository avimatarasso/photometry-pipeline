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
addpath('G:\code','G:\code\newPackage_020922')
% CUSTOMIZE


%% Only edit this section

% point to tanks (enter file path)
%path_to_data = 'G:\PhotomData\NEED TO ANALYZE\STIM ASSAYS\LCterminals\21060811_15Hz3sec_10Hz30sec\Avi_autoStim-210218-161112\dlight_CA1'; %must end in \
%savename = 'dlight_BLA';

path_to_data = 'G:\PhotomData\NEED TO ANALYZE\STIM ASSAYS\LCterminals\Dual GrabNE + dLight\21061416_10Hz3sec_15Hz30sec\Avi_autoStim-210218-161112\dlight_BLA\';
savename = 'dlight_BLA';
subjectPattern = ['*' savename '*'];
%subjectPattern = ['*31*'];

% Change subON to 1 if you want to subtract the LLS fit 
% Must customize if you would like to subtract 405 signal from 470 signal 
subON  = 1; %1 for both sensors
subFit = 1; %1 to subtract the linear least squares fit

% stimON and fearON should be 1 if you want the TTL from fear box or you
% want recorded stim points, leave as zero otherwise
% are you using stim or fearTTL? if you do both, you might need to adjust code
details.stimON = 1;
details.fearON = 0;

% do you want all figures to be shown? allfigures == 1 will show all
% check = 1 will show you raw figs, linear least squares fit, and df/f%
allfigures = 1;
check = 1;

%do not edit
saveON = 1; %save mat files
saveTXTON = 1; %save txt files 0  = won't save

% Bruchas lab specifics
details.static = 0; %1 if run on static system in Bruchas lab
details.auto   = 1; %1 if run with automated stim protocol - make sure tank name matches one below
%define tankname for tdt extraction
tankname = 'Avi_stim-210608-085426';% 'Avi_stim-210913-102052';
if details.auto
    tankname = 'Avi_autoStim-210907-110207'; % 'Avi_autoStim-210907-110207';
    %tankname = 'Avi_autoStim2-220504-084819';
    %tankname = 'Avi_autoStim_OneTrain-211104-151836';
    tankname = 'Avi_autoStim_OneTrain-220819-103307';
   % tankname = 'Avi_StimRecord-220209-101641';
    tankname = 'Avi_autoStim-210218-161112';
end

% What period of time do you want to baseline to? This will be your F0
BLlength = 160;

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
%workdir = workdir(3:end); %unmute this line if you don't want to use all
%files, CUSTOMIZE which

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

%convert all of the files 
if saveTXTON
    [dat470,dat405,Dts,timing] = tdtCONVERT(details, tankdir, tankname, blockname);
else
    [dat470,dat405,Dts] = tdtCONVERT(details, tankdir, tankname, blockname);
end
%sampling Frequency
Fs = round(1/(max(Dts)/length(Dts))); %sample frequency Hz

%% Clip the first second to get rid of artifact
dat405(1:Fs) = mean(dat405(1:BLlength*Fs));
dat470(1:Fs) = mean(dat470(1:BLlength*Fs));

%% preprocessing, fit to linear least squares/exponential, subtract isobestic
if subON
    fit405 = LLS(dat470,dat405);
%    fit405 = LLS(dat470(1:BLlength*Fs),dat405(1:BLlength*Fs));
    if subFit
        dataFilt = (dat470 - fit405)./fit405; %470 - fit(405) - F0 then normalized to 405fit
        dataFilt = dataFilt - mean(dataFilt(1:BLlength*Fs)); % Center your dF/F? CUSTOMIZE
    else
        dat470 = dat470 - dat405;  
        dataFilt = dat470/mean(dataFilt(1:BLlength*Fs));
        dataFilt = dataFilt - mean(dataFilt(1:BLlength*Fs)); % Center your dF/F? CUSTOMIZE
    end
    
    if check ==1
        % plot raw signals, linear least squares fit, and the df/f%
        if allfigures 
            figure
        else
            figure(1), clf
        end
        subplot(3,1,1)
        plot(Dts,dat470), hold on
        plot(Dts,dat405)
        legend('470','405');%'LLSfit405', 'df/f + mean(dat470)')
        ylim([220 280])
        title('raw traces')
        ylabel('mV')
        xlabel('time (s)')
        hold off

        subplot(3,1,2)
        plot(Dts,dat470), hold on
        plot(Dts,fit405)
        legend('470','llsfit(405)');%'LLSfit405', 'df/f + mean(dat470)')
        %ylim([mean(fit405)-30 mean(fit405)+30])
        title([filename(1:5) ' linear least squares fit - just LLS subtraction'])
        ylabel('mV')
        xlabel('time (s)')
        hold off
        
        subplot(3,1,3)
        plot(Dts,dataFilt*100)
        ylabel('dF/F %')
        xlabel('time (s)')
        ylim([-5 10])
        title('dF/F% of 470')
        hold off
    end
    
    % Do you want to subtract an exponential decay or LLS? 
    % Make sure the one you want is ending up in the dfF variable
    f2 = fit(Dts,dat470,'exp2');
    fitcurve= f2(Dts);
    dataFilt2 =  dat470 - fitcurve;
    
    dataFilt2 = (dat470 - dat405)./dat405; %470 - fit(405) - F0 then normalized to 405fit
    dataFilt2 = dataFilt2 - mean(dataFilt2(1:BLlength*Fs)); % Center your dF/F? CUSTOMIZE

    if check ==1
        % plot raw signals, linear least squares fit, and the df/f%
        if allfigures 
            figure
        else
            figure(1), clf
        end
        subplot(3,1,1)
        plot(Dts,dat470), hold on
        plot(Dts,dat405)
        legend('470','405');%'LLSfit405', 'df/f + mean(dat470)')
        ylim([220 280])
        title('raw traces')
        ylabel('mV')
        xlabel('time (s)')
        hold off

        subplot(3,1,2)
        plot(Dts,dat470), hold on
        plot(Dts,fit405)
        legend('470','470expfit');%'LLSfit405', 'df/f + mean(dat470)')
        ylim([mean(fit405)-30 mean(fit405)+30])
        title([filename(1:5) '470 exponential fit'])
        ylabel('mV')
        xlabel('time (s)')
        hold off
        
        subplot(3,1,3)
        plot(Dts,dataFilt2*100)
        ylabel('dF/F %')
        xlabel('time (s)')
        ylim([-5 10])
        title('dF/F% of (470 - 405)/405')
        hold off
    end
    
    
    
else
    %fit405 = LLS(dat470,dat405);
    %dat470 = dat470 - dat405;
    dataFilt = (dat470)./mean(dat470(1:BLlength*Fs));% - fit405)./fit405; %470 - fit(405) - F0 then normalized to 405fit
    dataFilt = dataFilt - mean(dataFilt(1:BLlength*Fs)); % Center your dF/F? CUSTOMIZE

    if check ==1
        % plot raw signals, linear least squares fit, and the df/f%
        if allfigures 
            figure
        else
            figure(1), clf
        end
        
        plot(Dts,dataFilt*100)
        ylabel('dF/F %')
        xlabel('time (s)')
        ylim([-5 10])
        title('dF/F% of 470, without 405 subtraction')
        hold off
    end
    
    
end
%% low pass filter

 %N = 1;  % order
 %cutoff_Hz = 10;  %
 %[b,a]=butter(N,cutoff_Hz/(Fs/2),'low');
 %dataFilt = filter(b,a,subdat);

%% calculate dF/F

dfF = dataFilt; %already dF/F %(dataFilt - median(dataFilt))./abs(median(dat470)); % this gives deltaF/F
dfFPerc = dfF.*100; % make a percentage
data1 = dfFPerc; %data1
        
% plot subtracted signal w/o baseline correction
%{
if allfigures 
    figure
else
    figure(6);
end
%plot(Dts,subdat,Dts,fitcurve);
plot(Dts,subdat);
title([filename(1:5) ': Subtracted signal'])
xlabel('time(s)')
ylabel('raw F')
%line([timingFear timingFear]', repmat([-30; 30],[1,length(timingFear)]),'Color','black')
%}
% plot final, baseline-corrected signal
%{
if allfigures 
    figure
else
    figure(7);
end
plot(Dts,data1);
title([filename(1:5) ': Signal with corrected baseline'])
xlabel('time(s)')
ylabel('deltaF/F')
ylim([-20 20])
%line([timingFear timingFear]', repmat([-30; 30],[1,length(timingFear)]),'Color','black')
%}

%% Recalculating the way Christian used to do

matlabFit = fit(Dts,(dat470-dat405),'exp2');
fitcurv = matlabFit(Dts);
dataFiltOld = dat470 - dat405 - fitcurv;
normDat1 = (dataFiltOld - median(dataFiltOld))./abs(median(dat470)); % this gives deltaF/F
normDat = normDat1.*100; % make a percentage
dataOld = normDat;
figure
plot(Dts,dataOld)



%% Save file as .mat file with specified filename
if saveON && saveTXTON
    if details.fearON && ~details.stimON
        save(filename,'dat470','dat405','fit405', 'dfF', 'dfFPerc','Dts','timing', 'f2', 'dataFilt2',"dataOld");
    elseif ~details.fearON && details.stimON
        save(filename,'dat470','dat405','fit405', 'dfF', 'dfFPerc','Dts','timing', 'f2', 'dataFilt2',"dataOld");
    elseif details.fearON && details.stimON
       save(filename,'dat470','dat405','fit405', 'dfF', 'dfFPerc','Dts','timingFear','timingStim', 'f2', 'dataFilt2',"dataOld");
    else
        save(filename,'dat470','dat405','fit405', 'dfF','dfFPerc','Dts', 'f2', 'dataFilt2',"dataOld");
    end  
elseif saveON
    if details.fearON && ~details.stimON
        save(filename,'dat470','dat405','fit405', 'dfF', 'dfFPerc','Dts', 'f2', 'dataFilt2',"dataOld");
    elseif ~details.fearON && details.stimON && subON
        save(filename,'dat470','dat405','fit405', 'dfF', 'dfFPerc','Dts', 'f2', 'dataFilt2',"dataOld");
    elseif ~details.fearON && details.stimON && ~subON
        save(filename,'dat470','dat405', 'dfF', 'dfFPerc','Dts');
    elseif details.fearON && details.stimON
       save(filename,'dat470','dat405','fit405', 'dfF', 'dfFPerc','Dts', 'f2', 'dataFilt2',"dataOld");
    else
        save(filename,'dat470','dat405','fit405', 'dfF','dfFPerc','Dts', 'f2', 'dataFilt2',"dataOld");
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
        
        %normal
        tmp2 = tmp{1}; tmp2 = strsplit(tmp2,'-'); 
        tmpDate = tmp{3}; tmpDate = strsplit(tmpDate,'-'); tmpDate = tmpDate{2}; 
        
        subStr = [tmp2{1} '-' tmp2{2}];
         region = savename;

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


