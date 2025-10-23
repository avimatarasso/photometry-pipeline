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