function finalSave(savename, sensor, region, allData,timingIdxs,timevec,path_to_data,dataForPlots)
% save all the variables and figures
 
if isfolder([region '_' sensor])
    cd([region '_' sensor])
    if exist(savename,'var')
        save(savename, 'allData','timingIdxs','timevec','path_to_data','dataForPlots','-append')
    else
        save(savename, 'allData','timingIdxs','timevec','path_to_data','dataForPlots')
    end
    cd('..')
else
    mkdir([region '_' sensor])
    cd([region '_' sensor])
    if exist(savename,'var')
        save(savename, 'allData','timingIdxs','timevec','path_to_data','dataForPlots','-append')
    else
        save(savename, 'allData','timingIdxs','timevec','path_to_data','dataForPlots')
    end
    cd('..')
    % close all
end

%% Save All Figures
answer = questdlg('All graphs are saved. Would you like to save the individual figures too?', ...
	'Saving?', ...
	'Yes','No','No');
switch answer
    case 'Yes'
        figureSaveName = savename;
        if exist([region '_' sensor],'dir')
            cd([region '_' sensor])
        else
            mkdir([region '_' sensor])
            cd([region '_' sensor])
        end
        tempdir = pwd;
        FolderName = tempdir;   % Your destination folder
        FigList = findobj(allchild(0), 'flat', 'Type', 'figure'); FigList2 = FigList(1:end);
        savefig(FigList2, fullfile(FolderName,[figureSaveName '.fig']));
       % save([savename '.mat']);
        disp('Figures have been saved!')        
        cd('..')
    case 'No'
        disp('You may manually save figures if you want.')
end

