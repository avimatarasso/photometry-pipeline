function [nextHeadEntry, numActual] = headEntryAdjust(rawtimings, eventLabel, subjectLabel, plot)

% this code creates a histogram for head entries for one photometry trace/
% experiment and looks at the histogram
% takes raw timings in column one of an array, and an event label in the
% second column
% subjectLabel can be taken out, i use it for titles. 
% plot == 1 will create a histogram for the subject you're looking at

eventTime = rawtimings(:,1);
eventtype = round(rawtimings(:,2));
eventTime = eventTime(eventtype==eventLabel); %head entry


if isempty(eventTime)
    warning([subjectLabel ' had no events!'])
    nextHeadEntry = [];
    numActual = 0;
else
    sucroseDel = rawtimings(:,1);
    eventtype = round(rawtimings(:,2));
    sucroseDel = sucroseDel(eventtype==1); %sucrose delivered here

    nextHeadEntry = zeros(length(sucroseDel),1);
    nbins = 300; % CUSTOMIZE 6 second bins because (30 minute(*60s) session)/nbins = 6 s bins

    idxsToDelete = [0; (diff(eventTime)<10)];
    eventTime = eventTime(~idxsToDelete);
    
    if plot && ~isempty(eventTime) 


        figure
        a = histogram(eventTime/60,nbins, 'DisplayStyle', 'stairs','EdgeColor','r');
        hold on
            xEnd = eventTime(end);
        ylabel('head entries')
        xlabel('time in minutes, 6 sec bins') %CUSTOMIZE
        maxVal = max(5*mod(max(a.Values),5),max(a.Values));
        ylim([0 (5-mod(maxVal,5))+maxVal])
        title(['distribution of head entries for ' subjectLabel]) %CUSTOMIZE

    end
    numActual = 0;
    for i = 1:length(sucroseDel)
        nextHeadEntries  = find((eventTime - sucroseDel(i) > 0)); %originally just first boolean
        %NEED TO GET RID OF THE LAST ONE IF IT IS THE SAME!    
        if ~isempty(nextHeadEntries) %need to make sure there was a head entry after the sucrose delivery!
            if i > 1 && numActual >=1
                lastHeadEntry = nextHeadEntry(i-1); %record the last head entry
            else 
                lastHeadEntry = 0; 
            end

            tempNextEntry = eventTime(nextHeadEntries(1));
            if tempNextEntry ~= lastHeadEntry || i == 1 % We don't want to accidentally record the same head entry twice
                nextHeadEntry(i) = tempNextEntry;
                numActual = numActual + 1; % 
            end
        end

        if plot
           L = line([sucroseDel(i) sucroseDel(i)]./60,([0 20]));
        end
    end
end
end