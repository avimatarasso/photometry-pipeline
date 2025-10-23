function createPlots(dataForPlots)
% You might want to customize each section of this code. It will create
% plots in the style that Avi typically likes them in

allData   = dataForPlots.allData;
indData   = dataForPlots.indData;
firstEvent = dataForPlots.firstEvent;
minTime   = dataForPlots.minTime;
timeBefore = dataForPlots.timeBefore;
timeAfter = dataForPlots.timeAfter;
dsFactor  = dataForPlots.dsFactor;
FS        = dataForPlots.FS;
titleName = dataForPlots.titleName;
workdir   = dataForPlots.workdir;
trials    = dataForPlots.trials;
savename  = dataForPlots.savename;
subName  = dataForPlots.subName;
trialsOI    = dataForPlots.trialsOI;
xlims = [-timeBefore timeAfter];
region  = dataForPlots.region;
sensor  = dataForPlots.sensor;

%% event triggered average
figure(61)
hold on

% make an x axis 
xlong = linspace(-timeBefore,timeAfter,length(allData));
x = decimate(linspace(-timeBefore,timeAfter,length(allData)),dsFactor);

% get rid of any nan points
if any(any(isnan(allData)))
    cols2getridof = isnan(allData(1,:));
    allData(:,cols2getridof) = [];
    warning([num2str(length(cols2getridof)) ' columns were deleted!'])
end
if size(allData, 2) > size(allData,1)
    allData = allData';
end

y = mean(allData,2); centerY = mean(mean(allData(xlong>-timeBefore & xlong<0,:)));
y = (y - centerY); %center the data
[peakY,ebIdx] = findpeaks(y'); [maxY,maxIdx] = max(peakY); ebIdx = ebIdx(maxIdx);
xMax = xlong(ebIdx); 
if y(ebIdx) ~= maxY
    warning('you may have different max value') 
end
y = decimate(double(y),dsFactor); %y = movmedian(y,6);
y = smoothdata(y,'SmoothingFactor',0.01);
sem = std(allData,0,2)./sqrt(size(allData,2)); % sem = std/sqrt(n)
semPeak = sem(ebIdx);
eb = decimate(double(sem'),dsFactor);
lineProps.col{1} = [0 0.5 0];
mseb(x,y,eb,lineProps,1);
%plot(linspace(-timewindow,timewindow,length(photoPerLick)),photoPerLick,'b')%,linspace(-timewindow,timewindow,length(photoPerLick)),LickTrig,'y')
%legend('SEM','Mean')
L = line([0 0],[-3 20]);
%L = line([30 30],[-3 5]);
if exist('xlims','var')
    xlim(xlims)
else
    xlim([-timeBefore timeAfter])
end
ylim([-1 5])
%axis([-timewindow timewindow -0.6 1.2])
set(L,'Color','black')
xlabel('Time aligned to event (s)')
ylabel('z-score')
title([titleName ', n = ' num2str(length(workdir))])
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
set(gca,'FontSize', 16)
xlabel('Time (s)', 'FontSize', 22)
ylabel('Z-score ','FontSize', 22)
%set(gcf, 'Position',  [100, 100, 600, 600])
set(gca,'FontName','Arial')
box off

hold off

if exist([region '_' sensor],'dir')
    cd([region '_' sensor])
else
    mkdir([region '_' sensor])
    cd([region '_' sensor])
end
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
saveas(gcf,[savename '_avged.svg'],'svg'), saveas(gcf,[savename '_avged.fig'],'fig'), 

save([savename '_peak_SEM.mat'],'maxY','semPeak','xMax')
cd('..')


%% WORK IN PROGRESS - early (all mice) to late (all mice) heatmap
figure(60)

% Until hold on, CUSTOMIZE for your dataset 
% i.e.: trials = 7    15     7     9     7    13     8    12     5
if size(trials,1)>1
    trials = trials';
end
cst = cumsum(trials); dcst = diff([1 cst]);
idxsToSort = round((timeBefore-timeBefore/4)*FS):round((timeAfter+timeBefore)/4*FS);

cstNew = [1 1+cst]; newIdxs = cstNew(1:end-1);
for ii = 1:max(dcst)-1
    newIdxs = [newIdxs cstNew+ii];
end
newIdxs = unique(newIdxs,'stable');
newIdxs(newIdxs>size(allData,2)) = [];
allChrono = allData(:,newIdxs);

hold on
% make a time vector for the heat map and plot the heat map
timevec=linspace(-timeBefore,timeAfter,length(allData));
imagesc(timevec,1:size(allChrono,2),allChrono')
L = line([0 0],[0 size(allChrono,1)+1]);
set(L,'Color','white')
set(L,'LineWidth',1)
xlabel('Time aligned to event (s)')
ylabel('Trial by mouse, bottom displays earliest trial')
%ytick(1:size(LickTrig,1))
title([titleName ' chronologically, n = ' num2str(length(workdir))])
cb = colorbar;
title(cb,'Z score')
caxis([-2 2])
colormap(flipud(brewermap([],'YlGnBu')))
colormap('parula')
if exist('xlims','var')
    xlim(xlims)
else
    xlim([-20 100])
end
ylim([0.5 size(allChrono,2)+0.5])
%set(gca, 'yticklabel', subName, 'ytick', cst);%trials/2 - 0.5 + (1:length(trialsOI):length(trialsOI)*length(subName)));
%for only one stim 
%set(gca, 'yticklabel', YLabelStrings, 'ytick', (1:length(YLabelStrings)));
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
set(gca,'FontSize', 16)
xlabel('Time (s)', 'FontSize', 22)
ylabel('Z-score ','FontSize', 22)
%set(gcf, 'Position',  [100, 100, 600, 600])
set(gca,'FontName','Arial')
box off

hold off

cd([region '_' sensor])

set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
saveas(gcf,[savename '_Chrono_heatmap.svg'],'svg'),saveas(gcf,[savename '_heatmap.fig'],'fig'),
cd('..')

%% ascend heatmap
figure(62)

%uncomment these if you prefer to get 
%[~,idx1] = sort(mean(LickTrig1(:,90000:150000),2),'ascend');
%LickTrig1 = LickTrig1(idx1,:);

cst = cumsum(trials);
idxsToSort = round((timeBefore-timeBefore/4)*FS):round((timeAfter+timeBefore)/4*FS);

%creating a matrix that has ascending peak response
[~,idx] = sort(mean(allData(idxsToSort,:),1),'ascend');
allData2=allData(:,idx);

hold on
% make a time vector for the heat map and plot the heat map
timevec=linspace(-timeBefore,timeAfter,length(allData));
imagesc(timevec,1:size(allData2,2),allData2')
L = line([0 0],[0 size(allData2,1)+1]);
set(L,'Color','white')
set(L,'LineWidth',1)
xlabel('Time aligned to event (s)')
ylabel('Trial by mouse, bottom displays earliest trial')
%ytick(1:size(LickTrig,1))
title([titleName ', n = ' num2str(length(workdir))])
cb = colorbar;
title(cb,'Z score')
caxis([-2 2])
colormap(flipud(brewermap([],'YlGnBu')))
colormap('parula')
if exist('xlims','var')
    xlim(xlims)
else
    xlim([-20 100])
end
ylim([0.5 size(allData2,2)+0.5])
%set(gca, 'yticklabel', subName, 'ytick', cst);%trials/2 - 0.5 + (1:length(trialsOI):length(trialsOI)*length(subName)));
%for only one stim 
%set(gca, 'yticklabel', YLabelStrings, 'ytick', (1:length(YLabelStrings)));
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
set(gca,'FontSize', 16)
xlabel('Time (s)', 'FontSize', 22)
ylabel('Z-score ','FontSize', 22)
%set(gcf, 'Position',  [100, 100, 600, 600])
set(gca,'FontName','Arial')
box off

load('cmapVar.mat')
colormap(BkRd)

hold off

cd([region '_' sensor])

set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
saveas(gcf,[savename '_ascendheatmap.svg'],'svg'),saveas(gcf,[savename '_heatmap.fig'],'fig'),

cd('..')

%% heat map (x dim: time, ydim: event, zdim: deltaF/F)
figure(65)

cst = cumsum(trials);
idxsToSort = round((timeBefore-timeBefore/4)*FS):round((timeAfter+timeBefore)/4*FS);
% uncomment if you want your data to be descending
%{
for ii = 1:length(workdir)

if ii == 1
    [~,idx] = sort(mean(allData(:,1:cst(ii)),1),'descend');
    allData2(:,1:cst(ii))=allData(:,idx);
else
    [~,idx] = sort(mean(allData(:,cst(ii-1)+1:cst(ii)),1),'descend');
    allData2(:,cst(ii-1)+1:cst(ii))=allData(:,cst(ii-1)+idx);
end13     5    16    15    18    11
end
%}
allData2 =allData;
hold on
timevec=linspace(-timeBefore,timeAfter,length(allData));
imagesc(timevec,1:size(allData2,2),allData2')
L = line([0 0],[0 size(allData2,1)+1]);
set(L,'Color','white')
set(L,'LineWidth',.25)

% plot horizontal lines
for ii = 1:length(workdir)
L2 = line([-150 size(allData2,1)],repmat(cst'+0.5, [1 2]));%([trials*ii + 0.5]));
set(L2,'Color','white')
set(L2,'LineWidth',.75)
end

xlabel('Time aligned to event (s)')
ylabel('Trial by mouse, bottom displays earliest trial')
%ytick(1:size(LickTrig,1))
title([titleName ', n = ' num2str(length(workdir))])
cb = colorbar;
title(cb,'Z score')
caxis([-2 2])
colormap(flipud(brewermap([],'YlGnBu')))
if exist('xlims','var')
    xlim(xlims)
else
    xlim([-20 100])
end
ylim([0.5 size(allData2,2)+0.5])
tickPts = cst;%[0 cst] + diff([1 cst 0])/2; tickPts= tickPts(1:end-1); 
set(gca, 'yticklabel', subName, 'ytick', tickPts); 
%for only one stim 
%set(gca, 'yticklabel', YLabelStrings, 'ytick', (1:length(YLabelStrings)));
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
set(gca,'FontSize', 16)
xlabel('Time (s)', 'FontSize', 22)
%ylabel('Subjects ','FontSize', 22)
%set(gcf, 'Position',  [100, 100, 600, 600])
set(gca,'FontName','Arial')
box off
load('cmapVar.mat')
colormap(BkRd)

hold off

cd([region '_' sensor])
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
saveas(gcf,[savename '_SubjectHeatmap.svg'],'svg'),saveas(gcf,[savename '_SubjectHeatmap.fig'],'fig'),

cd('..')

%% full session averaged, aligned to first event
figure(69)
hold on

% make an x axis 
xlong = 1/FS:1/FS:length(indData)/FS;

% get rid of any nan points
if any(any(isnan(indData)))
    cols2getridof = isnan(indData(1,:));
    indData(:,cols2getridof) = [];
    warning([num2str(length(cols2getridof)) ' columns were deleted!'])
end
if size(indData, 2) > size(indData,1)
    indData = indData';
end

adjustSamps = floor((FS*(firstEvent - min(firstEvent))));
tmpArr = [];
for iii=1:size(subName,2)
    addedSamp = adjustSamps(iii);
    if (floor(minTime*FS)+addedSamp)>length(indData)
        tmp = indData(addedSamp+1:floor(minTime*FS),iii);
    else
        tmp = indData(addedSamp+1:floor(minTime*FS)+addedSamp,iii);
    end
    if iii == 1
        tmpArr = [tmpArr tmp];
    else
        if length(tmpArr)>length(tmp) && ~isempty(tmp)
                tmpArr = [tmpArr(1:length(tmp),:) tmp];
            elseif length(tmp)>length(tmpArr)
                tmpArr = [tmpArr tmp(1:length(tmpArr))];
            else
                tmpArr = [tmpArr tmp];
        end
        if isempty(tmp)
            subName{iii} = [];
        end
    end
end
nData = tmpArr;
startingSamp = floor(min(firstEvent)*FS);
tmpx = xlong(floor(startingSamp:floor(minTime*FS)));
x = decimate(tmpx,dsFactor);

y = mean(nData,2); centerY = mean(mean(nData(startingSamp:round( FS*min(firstEvent)),:)));
ylong = (y - centerY); %center the data
[peakY,ebIdx] = findpeaks(y'); [maxY,maxIdx] = max(peakY); ebIdx = ebIdx(maxIdx);
xMax = xlong(ebIdx); 
if y(ebIdx) ~= maxY
    warning('you may have different max value') 
end

y = decimate(double(ylong(floor(startingSamp:end))),dsFactor); %y = movmedian(y,6);
y = smoothdata(y,'SmoothingFactor',0.01);
sem = std(nData,0,2)./sqrt(size(nData,2)); % sem = std/sqrt(n)
semPeak = sem(ebIdx);
if length(sem)>length(ylong)
    sem = sem(floor(min(firstEvent)*FS):floor(minTime*FS));
end
eb = decimate(double(sem'),dsFactor);
lineProps.col{1} = [0 0.5 0];
mseb(x(1:length(y)),y,eb(1:length(y)),lineProps,1);
%plot(linspace(-timewindow,timewindow,length(photoPerLick)),photoPerLick,'b')%,linspace(-timewindow,timewindow,length(photoPerLick)),LickTrig,'y')
%legend('SEM','Mean')
L = line([0 0],[-3 20]);
%L = line([30 30],[-3 5]);
xlimEnd = floor(length(nData)/FS);
xlim([startingSamp/FS xlimEnd])
ylim([-1 3])
%axis([-timewindow timewindow -0.6 1.2])
set(L,'Color','black')
xlabel('Time aligned to event (s)')
ylabel('z-score')
title([titleName ' full session aligned to first event, n = ' num2str(length(workdir))])
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
set(gca,'FontSize', 16)
xlabel('Time (s)', 'FontSize', 22)
ylabel('Z-score ','FontSize', 22)
%set(gcf, 'Position',  [100, 100, 600, 600])
set(gca,'FontName','Arial')
box off

hold off

if exist([region '_' sensor],'dir')
    cd([region '_' sensor])
    if exist('fullSess','dir')
        cd('fullSess')
    else
        mkdir('fullSess')
        cd('fullSess')
    end
else
    mkdir([region '_' sensor])
    cd([region '_' sensor])
end
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
saveas(gcf,[savename '_avgedSession1stevnt.svg'],'svg'), saveas(gcf,[savename '_avgedSession1stevnt.fig'],'fig'), 
cd('../..')


figure(34)
cst = cumsum(trials);
idxsToSort = round((timeBefore-timeBefore/4)*FS):round((timeAfter+timeBefore)/4*FS);

%creating a matrix that has ascending peak response
[~,idx] = sort(mean(indData,1),'ascend');

hold on
% make a time vector for the heat map and plot the heat map
imagesc(xlong,1:size(indData,2),indData(:,idx)')
L = line([0 0],[0 size(indData,1)+1]);
set(L,'Color','white')
set(L,'LineWidth',1)
xlabel('Time aligned to event (s)')
ylabel('Trial by mouse, bottom displays earliest trial')
%ytick(1:size(LickTrig,1))
title([titleName ' fullSession, n = ' num2str(length(workdir))])
cb = colorbar;
title(cb,'Z score')
caxis([-2 2])

colormap(flipud(brewermap([],'YlGnBu')))
colormap('parula')
xlim([startingSamp/FS xlimEnd])
ylim([0.5 size(indData,2)+0.5])
%set(gca, 'yticklabel', subName, 'ytick', cst);%trials/2 - 0.5 + (1:length(trialsOI):length(trialsOI)*length(subName)));
%for only one stim 
%set(gca, 'yticklabel', YLabelStrings, 'ytick', (1:length(YLabelStrings)));
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
set(gca,'FontSize', 16)
xlabel('Time (s)', 'FontSize', 22)
ylabel('Z-score ','FontSize', 22)
%set(gcf, 'Position',  [100, 100, 600, 600])
set(gca,'FontName','Arial')
box off

hold off

cd([region '_' sensor])
cd('fullSess')
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
saveas(gcf,[savename '_fullSess1stevntheatmap.svg'],'svg'),saveas(gcf,[savename '_fullSess1stevntheatmap.fig'],'fig'),
cd('../..')

%% full session averaged, aligned to full session 
figure(70)
hold on

% make an x axis 
xlong = 1/FS:1/FS:length(indData)/FS;

% get rid of any nan points
if any(any(isnan(indData)))
    cols2getridof = isnan(indData(1,:));
    indData(:,cols2getridof) = [];
    warning([num2str(length(cols2getridof)) ' columns were deleted!'])
end
if size(indData, 2) > size(indData,1)
    indData = indData';
end

%adjustSamps = floor((FS*(firstEvent - min(firstEvent))));
tmpArr = [];
for iii=1:size(indData,2)
    tmp = indData(1:floor(minTime*FS),iii);
    if iii == 1
        tmpArr = [tmpArr tmp];
    else
        if length(tmpArr)>length(tmp) 
                tmpArr = [tmpArr(1:length(tmp),:) tmp];
            elseif length(tmp)>length(tmpArr)
                tmpArr = [tmpArr tmp(1:length(tmpArr))];
            else
                tmpArr = [tmpArr tmp];
        end
    end
end
nData = tmpArr;

tmpx = xlong(1:floor(minTime*FS));
x = decimate(tmpx,dsFactor);

y = mean(nData,2); centerY = mean(mean(nData(1:round( FS*min(firstEvent)),:)));
ylong = (y - centerY); %center the data
[peakY,ebIdx] = findpeaks(y'); [maxY,maxIdx] = max(peakY); ebIdx = ebIdx(maxIdx);
xMax = xlong(ebIdx); 
if y(ebIdx) ~= maxY
    warning('you may have different max value') 
end

y = decimate(double(ylong(1:end)),dsFactor); %y = movmedian(y,6);
y = smoothdata(y,'SmoothingFactor',0.01);
sem = std(nData,0,2)./sqrt(size(nData,2)); % sem = std/sqrt(n)
semPeak = sem(ebIdx);
if length(sem)>length(ylong)
    sem = sem(1:floor(minTime*FS));
end
eb = decimate(double(sem'),dsFactor);
lineProps.col{1} = [0 0.5 0];
mseb(x(1:length(y)),y,eb(1:length(y)),lineProps,1);
%plot(linspace(-timewindow,timewindow,length(photoPerLick)),photoPerLick,'b')%,linspace(-timewindow,timewindow,length(photoPerLick)),LickTrig,'y')
%legend('SEM','Mean')
L = line([0 0],[-3 20]);
%L = line([30 30],[-3 5]);
xlimEnd = floor(length(nData)/FS);
xlim([1 xlimEnd])
ylim([-1 3])
%axis([-timewindow timewindow -0.6 1.2])
set(L,'Color','black')
xlabel('Time aligned to event (s)')
ylabel('z-score')
title([titleName ' full session completely aligned, n = ' num2str(length(workdir))])
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
set(gca,'FontSize', 16)
xlabel('Time (s)', 'FontSize', 22)
ylabel('Z-score ','FontSize', 22)
%set(gcf, 'Position',  [100, 100, 600, 600])
set(gca,'FontName','Arial')
box off

hold off

if exist([region '_' sensor],'dir')
    cd([region '_' sensor])
    if exist('fullSess','dir')
        cd('fullSess')
    else
        mkdir('fullSess')
        cd('fullSess')
    end
else
    mkdir([region '_' sensor])
    cd([region '_' sensor])
end
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
saveas(gcf,[savename '_avgedSession.svg'],'svg'), saveas(gcf,[savename '_avgedSession.fig'],'fig'), 
cd('../..')


figure(35)
cst = cumsum(trials);
idxsToSort = round((timeBefore-timeBefore/4)*FS):round((timeAfter+timeBefore)/4*FS);

%creating a matrix that has ascending peak response
[~,idx] = sort(mean(indData,1),'ascend');

hold on
% make a time vector for the heat map and plot the heat map
imagesc(xlong,1:size(indData,2),indData(:,idx)')
L = line([0 0],[0 size(indData,1)+1]);
set(L,'Color','white')
set(L,'LineWidth',1)
xlabel('Time aligned to event (s)')
ylabel('Trial by mouse, bottom displays earliest trial')
%ytick(1:size(LickTrig,1))
title([titleName ' fullSession aligned in full, n = ' num2str(length(workdir))])
cb = colorbar;
title(cb,'Z score')
caxis([-2 2])

colormap(flipud(brewermap([],'YlGnBu')))
colormap('parula')
xlim([min(firstEvent) xlimEnd])
ylim([0.5 size(indData,2)+0.5])
%set(gca, 'yticklabel', subName, 'ytick', cst);%trials/2 - 0.5 + (1:length(trialsOI):length(trialsOI)*length(subName)));
%for only one stim 
%set(gca, 'yticklabel', YLabelStrings, 'ytick', (1:length(YLabelStrings)));
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
set(gca,'FontSize', 16)
xlabel('Time (s)', 'FontSize', 22)
ylabel('Z-score ','FontSize', 22)
%set(gcf, 'Position',  [100, 100, 600, 600])
set(gca,'FontName','Arial')
box off

hold off

cd([region '_' sensor])
cd('fullSess')
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
saveas(gcf,[savename '_fullSessheatmap.svg'],'svg'),saveas(gcf,[savename '_fullSessheatmap.fig'],'fig'),
cd('../..')

%% plot only some events, such as the last stims

%{
% event triggered average
figure(63)
allData = ad(:,[2 3 5 6 8 9]);
hold on
xlong = linspace(-timeBefore,timeAfter,length(nData));
x = decimate(linspace(-timeBefore,timeAfter,length(nData)),dsFactor);
if any(any(isnan(allData)))
    cols2getridof = isnan(allData(1,:));
    allData(:,cols2getridof) = [];
    warning('some columns were deleted!')
end
y = mean(allData,2); centerY = mean(mean(allData(1:round(timeBefore*FS),:)));
y = (y - centerY); %center the data
[peakY,ebIdx] = findpeaks(y'); [maxY,maxIdx] = max(peakY); ebIdx = ebIdx(maxIdx);
xMax = xlong(ebIdx); 
if y(ebIdx) ~= maxY
    warning('you may have different max value') 
end
y = decimate(y,dsFactor);
%y = smoothdata(y,'SmoothingFactor',0.1);
sem = std(allData,0,2)./sqrt(size(allData,2)); % sem = std/sqrt(n)
semPeak = sem(ebIdx);
eb = decimate(sem',dsFactor);
lineProps.col{1} = [0 0.5 0];
mseb(x,y,eb,lineProps,1);
%plot(linspace(-timewindow,timewindow,length(photoPerLick)),photoPerLick,'b')%,linspace(-timewindow,timewindow,length(photoPerLick)),LickTrig,'y')
%legend('SEM','Mean')
L = line([0 0],[-3 20]);
%L = line([30 30],[-3 5]);
if exist('xlims','var')
    xlim(xlims)
else
    xlim([-20 timeAfter])
end
ylim([-1 5])
%axis([-timewindow timewindow -0.6 1.2])
set(L,'Color','black')
xlabel('Time aligned to event (s)')
ylabel('z-score')
title([titleName ', n = ' num2str(length(workdir))])
hold off

if exist([region '_' sensor],'dir')
    cd([region '_' sensor])
else
    mkdir([region '_' sensor])
    cd([region '_' sensor])
end
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
saveas(gcf,[savename '_avged_lastStims.svg'],'svg'), saveas(gcf,[savename '_avged_lastStims.fig'],'fig'), 
saveas(gcf,[savename '_avged_lastStims.svg'],'svg'), saveas(gcf,[savename '_avged_lastStims.fig'],'fig'), 
saveas(gcf,[savename '_avged_lastStims.svg'],'svg'), saveas(gcf,[savename '_avged_lastStims.fig'],'fig'), 

save([savename '_peak_lastStims_SEM.mat'],'maxY','semPeak','xMax')
cd('..')
%}


end