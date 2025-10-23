function [nData] = centAndNormData(nData, FS, timeBefore)
%loop through timings

%find BLStart
BLend = ceil(timeBefore*FS);

%center around baseline
nData = nData - mean(nData(1:BLend,:),1);

%{
%endIdx   = startIdx + lengthOfBL;
bl  = nData(1:BLendIdx,:); %find baseline 
meanBL   = abs(mean(bl,1));
stdBL    = abs(std(bl,1));

%normalize the data after alignTime from the start of timeBefore
cData = nData-meanBL;
ncData= nData/stdBL;
%}
end