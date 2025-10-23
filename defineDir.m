function workdir = defineDir(path_to_data, commonStr)

for pp = 1:length(path_to_data)
% Organize the path, the files, and get rid of spaces for trialname
tempPath = path_to_data{pp};
cd([tempPath]), cd('..')
filesep = '\'; loc = strfind(tempPath, filesep); loc = loc(end-1:end);
folderName   = tempPath(loc(1)+1:loc(2)-1);
spaceLocs = strfind(folderName, ' ');  
%trialname = tempPath(loc(1)+1:loc(2)-1); trialname(spaceLocs) = '_';

cd(tempPath)

if ~exist('workdir','var')
    workdir = dir(['*' commonStr '*.mat']);%dir(['*' region '*.mat']);
else
    wdInitialLength = length(workdir);
    tempWD  = dir(['*' commonStr '*.mat']);
    tempWDlen = length(tempWD);
    workdir(wdInitialLength+1:tempWDlen+wdInitialLength) = tempWD; 
end
    
if isempty(workdir)
    workdir = dir('*.mat'); %[subStr '*.mat']);
    if isempty(workdir)
        warning('you have to change your workdir variable')
    end
end
end


end