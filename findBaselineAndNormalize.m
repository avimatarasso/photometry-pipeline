function [nData, bl, timingIdxs] = findBaselineAndNormalize(data, FS, timings, timeBefore, timeAfter)
% find baseline right before event and normalize to that
% timings [#timePts x 1] should just be the timings in seconds

lengthOfBL = floor(abs(timeBefore)*FS);
lengthOfData = floor((timeAfter+timeBefore)*FS);

%loop through htimings
for i = 1:length(timings)
    startIdx = ceil((timings(i)-timeBefore)*FS);
    endIdx   = startIdx + lengthOfBL;
    if endIdx > length(data)
        warning('Your timeAfter variable goes past the length of the data.')
    end
    bl(:,i)  = data(startIdx:endIdx-1); %find baseline 
    meanBL   = abs(mean(bl(:,i)));
    stdBL    = abs(std(bl(:,i)));

    %normalize the data after alignTime from the start of timeBefore
    endIdx   = startIdx + lengthOfData;
    nData(:,i) = (data(startIdx:endIdx-1)-meanBL)/stdBL;
    %nData(:,i) = data(startIdx:endIdx-1);
    timingIdxs(i,:)  = [startIdx endIdx]; 
end

end %end func

