function [subStr, subNumb, subjectLabel, txtName] = defineSubject(photoname, region,subjectIdx)
% define names of subject and text files
% THIS FUNCTION IS VERY INDIVIDUALIZED

    firstDash = findstr(photoname,'-'); firstDash = firstDash(1);
    subStr = [photoname(1:firstDash-1) '-']; 
    subNumb   = photoname(strfind(photoname, subStr) + length(subStr)); % will get sub
    subjectLabel = [subStr subNumb]; 
    
    %customize special cases with if statements
    if isempty(str2num(subStr(end))) %isempty(str2num(subNumb))
        idx1 = strfind(photoname, '-') + 1; idx1=idx1(1);
        subNumb   = photoname(idx1); % will get subject number  
        if strcmp(photoname(3),'-')
            subStr = [photoname(1:2) '-'];
        end
        subjectLabel = photoname(1:idx1); 
        %subNumb   = photoname(strfind(photoname, subStr) + length(subStr))-1; % will get sub
    end
    
    %check
    dateLabel = strsplit(photoname,'-'); dateLabel=dateLabel{3};

    txtName = dir([subjectLabel '*' region '*' dateLabel '*.txt']); 

    %%% Make sure the Timings you use are consistent, or change them every time
    % CUSTOMIZE 
    
    if isempty(txtName)
        subName{subjectIdx} = subjectLabel;
        txtName    = [subStr 'Stim_timesM' subNumb '.txt'];    
    
        %txtName = [subStr subjectLabel(end) '_' region '.txt'];
        txtName = [subjectLabel '_' region '.txt']; %txtName(7) = '_';
        try dir([txtName(1:end-4) '*.txt'])
            txtName = dir([txtName(1:end-4) '*.txt']); txtName=txtName(1).name;
        end
    end
    try txtName = txtName.name;
    end
    if ~isfile(txtName)
    %ELENA NAMIGN CONVENTION
        try txtName = [subName{subjectIdx} '*.txt'];
            txtName = dir([txtName(1:end-4) '*.txt']); txtName=txtName(1).name;
        end
    end
    
    if ~isfile(txtName)
        %txtName = [subStr subjectLabel(end) '_' region '.txt'];
        txtName = [subjectLabel '_' region '*.txt'];  
        txtName = dir(txtName); 
        if ~isempty(txtName)
            if isfile(txtName)
            txtName = txtName(1).name;       
            end
        end
    end

    if ~isfile(txtName)
        delims = regexp(photoname,'-'); delim2 = delims(2);
       tmptxt = [photoname(1:delim2-1) '*.txt']; tmptxt(delims(1))='_';

        txtName = dir(tmptxt);
        txtName = txtName(1).name;
    end
    
    photoname
    txtName
end