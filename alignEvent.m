function [nData, timingIdxs, sampdiff] = alignEvent(data, FS, timings, timeBefore, timeAfter)
% find baseline of each time point
% timings [#timePts x 1] should just be the timings in seconds

lengthOfData = floor((timeAfter+timeBefore)*FS);
nData = nan(lengthOfData,length(timings));
timingIdxs = nan(2,length(timings))';
%loop through timings
for i = 1:length(timings)
    startIdx = ceil((timings(i)-timeBefore)*FS);
    if startIdx == 0 
        startIdx = 1;
    end

    endIdx   = startIdx + lengthOfData;
    %Align data based on time before and after event
    if startIdx>=1 && endIdx<length(data)
        nData(:,i) = data(startIdx:endIdx-1);
        timingIdxs(i,:)  = [startIdx endIdx]; 
    else
        continue
    end
end
sampdiff = timings*FS - timingIdxs(:,1); %sample difference from start of BL
sampdiff = sampdiff(1);
timeAx = -timeBefore:1/FS:timeAfter-1/FS;
    
end 
