% Code written by Christian Pedersen
% Michael Bruchas Lab - UW

% Extracts photometry data from TDT Data Tank and outputs time series as
% .mat file

% REQUIRES tdt2mat.m file to be in same folder (filepath)

%% Reset MatLab workspace - clears all variables and the command window

clear all;  % clear all variables
close all;  % close all open graphs
clc   % clear command window

addpath('F:\code') %D:\MATLAB\Eric scripts %change this


%% Only edit this section
      
% point to tanks (enter file path)
path_to_data = 'F:\PhotomData\NEED TO ANALYZE\OPERANT\photom\'; %dont include tank %change this
subjectPattern = '1*';

%path_to_data = 'G:\PhotomData\dLight\cocaine+stim_111820\'; %dont include block
%subjectPattern = '139*';

% Change subON to 1 if you want to subtract 405 signal from 470 signal 
% stimON and fearON should be 1 if you want the TTL from fear box or you
% want record
%leave as zero otherwise
subON  = 1; %1 for both sensors
stimON = 0;
saveON = 1; %save mat files
saveTXTON = 0; %save txt files 0  = won't save

%%
%cd(path_to_data)
% set up Tank name, variables to extract//
tankdir = [path_to_data];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comment out which one you aren't looking at
tankname = 'Raaj_twomice-210204-181611';
tankname = 'Daniel_2mice-190722-131320'; 
%tankname = 'Daniel_2mice-210910';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

workPath = [path_to_data '\' tankname];
cd(workPath)
workdir = dir(subjectPattern); %dir(pwd) for all directories
%workdir = workdir(3:end); %unmute this line if you don't want to use
%pattern
% Get a logical vector that tells which is a directory.
dirFlags = [workdir.isdir];
% Extract only those that are directories.
workdir = workdir(dirFlags);

for ii = 1:2 %9:length(workdir)
%tdt files must match protocol
blockname = workdir(ii).name; % name of your file
if ~isempty(strfind(blockname,'.mat'))
    blockname = blockname(1:end-4); %get rid of .mat at the end
end
% pick any file name to save time series (must end in .mat)
filename = [blockname '.mat'];

%for debugging when you already have some files
%if exist(filename)==2
%    continue
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if contains(tankname,'Daniel')
    storenames = {'470A'}; % name of stores to extract from TDT (usu. 4-letter code) 
    %LMag is the demodulated data, may also have other timestamps etc

    storenames_iso = {'405A'};

    storenames2     = {'470B'}; %stim

    storenames2_iso = {'405B'}; % Fear 
else
    storenames = {'465A'}; % name of stores to extract from TDT (usu. 4-letter code) 
    %LMag is the demodulated data, may also have other timestamps etc

    storenames_iso = {'405A'};

    storenames2     = {'465C'}; %stim

    storenames2_iso = {'405C'}; % Fear 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract

for k = 1:numel(storenames)
  storename = storenames{k};
  S{k} = tdt2mat(tankdir, tankname, blockname, storename);
  
  storenames_iso = storenames_iso{k};
  S_iso{k} = tdt2mat(tankdir, tankname, blockname, storenames_iso);

  storename2 = storenames2{k};
  S2{k} = tdt2mat(tankdir, tankname, blockname, storename2);
  
  storename2_iso = storenames2_iso{k};
  S2_iso{k} = tdt2mat(tankdir, tankname, blockname, storename2_iso);
end


% Massage data and get time stamps

LMag = S{1}; %add more if you extracted more stores above
% LMag2 = S{2};
% For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
chani1 = LMag.channels==1;
chani2 = LMag.channels==2;

LMag2 = S_iso{1}; %add more if you extracted more stores above
% LMag2 = S{2};
% For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
chani21 = LMag2.channels==1;
chani22 = LMag2.channels==2;

LMag3 = S2{1}; %add more if you extracted more stores above
% LMag3 = S3{2};
% For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
chani31 = LMag3.channels==1;
chani32 = LMag3.channels==2;

LMag4 = S2_iso{1}; %add more if you extracted more stores above
% LMag3 = S3{2};
% For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
chani41 = LMag4.channels==1;
chani42 = LMag4.channels==2;

% Get LMag data as a vector (repeat for each channel)
dat470A = LMag.data(chani1,:);
dat470A = reshape(dat470A', [],1); % unwrap data from m x 256 array
% dat405A = LMag.data(chani21,:);
% dat405A = reshape(dat405A', [],1); % unwrap data from m x 256 array

% Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
ts = LMag.timestamps(chani1);
t_rec_start = ts(1);
dt = datetime( t_rec_start, 'ConvertFrom', 'posixtime' );

ts = ts-ts(1); % convert from Unix time to 'seconds from block start'
ts = bsxfun(@plus, ts(:), (0:LMag.npoints-1)*(1./LMag.sampling_rate));
ts = reshape(ts',[],1);

%%%%%%%%%%%%%%%%%%
 
dat405A = LMag2.data(chani21,:);
dat405A = reshape(dat405A', [],1); % unwrap data from m x 256 array
% dat405A = LMag.data(chani21,:);
% dat405A = reshape(dat405A', [],1); % unwrap data from m x 256 array

% Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
ts2 = LMag2.timestamps(chani21);
t_rec_start2 = ts2(1);

ts2 = ts2-ts2(1); % convert from Unix time to 'seconds from block start'
ts2 = bsxfun(@plus, ts2(:), (0:LMag2.npoints-1)*(1./LMag2.sampling_rate));
ts2 = reshape(ts2',[],1);

% 465 C and 405 C below
dat465C = LMag3.data(chani31,:);
dat465C = reshape(dat465C', [],1); % unwrap data from m x 256 array

% Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
ts3 = LMag3.timestamps(chani31);
t_rec_start = ts3(1);
dt = datetime( t_rec_start, 'ConvertFrom', 'posixtime' );
ts3 = bsxfun(@plus, ts3(:), (0:LMag.npoints-1)*(1./LMag.sampling_rate));
ts3 = reshape(ts3',[],1);

dat405C = LMag4.data(chani41,:);
dat405C = reshape(dat405C', [],1); % unwrap data from m x 256 array

ts4 = LMag4.timestamps(chani41);
t_rec_start4 = ts4(1);
dt4 = datetime( t_rec_start4, 'ConvertFrom', 'posixtime' ); 

ts4 = ts4-ts4(1); % convert from Unix time to 'seconds from block start'
ts4 = bsxfun(@plus, ts4(:), (0:LMag4.npoints-1)*(1./LMag4.sampling_rate));
ts4 = reshape(ts4',[],1);
timingStart = seconds(dt4 - dt); 



%% Smooth signal

dat470A = smoothdata(dat470A,'SmoothingFactor',0.1);
dat405A = smoothdata(dat405A,'SmoothingFactor',0.1);
dat465C = smoothdata(dat465C,'SmoothingFactor',0.1);
dat405C = smoothdata(dat405C,'SmoothingFactor',0.1);
    
a = [size(dat470A,1),size(dat405A,1)]; [minA, minIdx] = min(a);
c = [size(dat465C,1),size(dat405C,1)]; [minC, minIdx] = min(c);
dat470A = dat470A(1:minA); dat405A = dat405A(1:minA);
dat470C = dat465C(1:minC); dat405C = dat405C(1:minC);
ts = ts(1:minA); ts2 = ts2(1:minA);
ts3 = ts3(1:minC); ts4 = ts4(1:minC);

%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subtracted 405 signal from 470 signal % raw final signal
subdat = dat470A; %-dat405A;

if subON 
    subdat = dat470A-dat405A;
    subdat2= dat465C-dat405C;
end

Dts=ts;
Dts2=ts3;

% plot raw 405 ch, raw 470 ch and subtracted signal
figure(2), clf
hold on
plot(ts2,dat405A,'r');
plot(ts,dat470A,'b');
title('both 405 and 470 signals')
%plot(ts,subdat,'g');
xlabel('time(s)')
ylabel('amplitude')
a=find(Dts>5);
dataTot = [dat470A(a:end); dat405A(a:end)];
maxdT=max(dataTot);
mindT=min(dataTot);
ylim([mindT maxdT])
legend('405','470')
hold off


% Bandpass filter data 
Fs = round(1/(max(Dts)/length(Dts))); %sample frequency Hz

%% IF PHOTOMETRY RECORDING BASELINE DROPS OFF:

%7200 cutoff P6-3 PR 1
% subdat = subdat(1:Fs*7200);
% Dts = Dts(1:length(subdat));

% % 2800 cutoff P11-1 PR 2
% subdat = subdat(1:Fs*2850);
% Dts = Dts(1:length(subdat));

%% instead of HPF: fit 2nd order exponential and subtract

if subON
    %subdat = subdat - mean(subdat);
    f2 = fit(Dts,dat405A,'poly2');
    fit405A= f2(Dts);
    dataFilt =  dat470A - fit405A;

    subdat2 = subdat2 - mean(subdat2);
    f2 = fit(Dts2,dat405C,'poly2');
    fit405C= f2(Dts2);
    dataFilt2 =  dat465C - fit405C;
else
    dat470A = dat470A - mean(dat470A);
    f2 = fit(Dts,dat470A,'exp2');
    fitcurve= f2(Dts);
    dataFilt =  dat470A - fitcurve;
    
    dat465C = dat465C - mean(dat465C);
    f2 = fit(Dts,dat465C,'exp2');
    fitcurve= f2(Dts);
    dataFilt2 =  dat465C - fitcurve;
end

%% low pass filter

 %N = 1;  % order
 %cutoff_Hz = 10;  %
 %[b,a]=butter(N,cutoff_Hz/(Fs/2),'low');
 
 %dataFilt = filter(b,a,subdat);

%% calculate dF/F

%
normDat1 = (dataFilt - median(dataFilt))./abs(median(dat470A)); % this gives deltaF/F
normDat = normDat1.*100; % make a percentage

dat470A = normDat;

normDat1 = (dataFilt2 - median(dataFilt2))./abs(median(dat465C)); % this gives deltaF/F
normDat = normDat1.*100; % make a percentage

dat465C = normDat;

if subON
    normDat1 = (dataFilt - median(dataFilt))./abs(median(subdat)); % this gives deltaF/F
    normDat = normDat1.*100; % make a percentage

    subdat = normDat;

    normDat1 = (dataFilt2 - median(dataFilt2))./abs(median(subdat2)); % this gives deltaF/F
    normDat = normDat1.*100; % make a percentage

    subdat2 = normDat;
end
%

% plot subtracted signal w/o baseline correction
if exist('fitcurve')
figure(6);
    plot(Dts,subdat,Dts,fitcurve(1:length(Dts)));
title('Subtracted signal')
xlabel('time(s)')
ylabel('raw F')
end

% plot final, baseline-corrected signal
figure(7);
plot(Dts,dat470A);
title('Signal with corrected baseline')
xlabel('time(s)')
ylabel('deltaF/F')
ylim([-20 20])

%%

A465 = dat470A;
A405 = dat405A; 
C465 = dat465C;
C405 = dat405C; 
if subON
   AsubIso = subdat;  %box 6 for Raaj, box 1 for daniel
   CsubIso = subdat2; %box 5 for Raaj, box 2 for daniel
end

%% Save file as .mat file with specified filename
if saveON
if subON
    save(filename,'AsubIso','CsubIso','dt','Dts');
else
    save(filename,'A465','C465','dt','Dts');
end
end

end


%% RUN THIS SECTION SEPARATELY!

%naming conventions for data 
% If in daniel, first is box 1, second is box 2
% If in raaj, first is box 5 (BOX D - FIBPHO2), second is box 6 (BOX C - FIBPHO1)

workDir = 'F:\PhotomData\NEED TO ANALYZE\OPERANT\FROM WEEK 2\BLA,CA1 operant\Daniel_2mice-210910'; %path_to_data;
%workDir = 'F:\PhotomData\NEED TO ANALYZE\MAGTRAIN042021\Raaj_twomice-210204-181611';
cd(workDir)
tmp = strsplit(workDir,'\'); tmpDir = tmp{end};
    allFiles = dir('*.mat');

for i = 1:length(allFiles)
    tmpName = allFiles(i).name;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% CUSTOMIZE
    deconstructed = regexp(tmpName,'\d*','Match');
    decon = deconstructed(1:end-2);
    deconDay = deconstructed{end-1};
    if length(decon) == 3
        Str1 = [decon{1} '-' decon{2} '_' deconDay '.mat'];
        Str2 = [decon{1} '-' decon{3} '_' deconDay '.mat'];
    elseif length(decon) == 4
        Str1 = [decon{1} '-' decon{2} '_' deconDay '.mat'];
        Str2 = [decon{3} '-' decon{4} '_' deconDay '.mat'];    
        Str1 = [decon{1} '-' decon{2} '_' deconDay '.mat'];
        Str2 = [decon{1} '-' decon{3} '_' deconDay '.mat'];    
    else
        warning('You may want to check your files match the naming convention set')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    load(tmpName)
    %for choosing whether A or C is which text file
    if contains(tmpDir,'Daniel')
        data1 = AsubIso;
        data2TMP = CsubIso;
    elseif contains(tmpDir,'Raaj')
        data2TMP = AsubIso;
        data1 = CsubIso;        
    else
        warning('are you in the right location?')
    end
    
    if ~exist('separatedFiles','dir')
        mkdir('separatedFiles')
        cd('separatedFiles')
    else
        cd('separatedFiles')
    end
    
    save(Str1,'data1','dt','Dts');
    data1 = data2TMP; %make the naming convention the same for all files 
    save(Str2,'data1','dt','Dts');
    
    cd(workDir)
    clearvars -except tmpDir workDir allFiles path_to_data
end
    
