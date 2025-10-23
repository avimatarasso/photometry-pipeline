function plotIndividual(figureN, eventTime, raw, data1, subName, subjectIdx, dsFactor, time)
    
    raw405 = raw.rawIso(subjectIdx,:);
    raw470 = raw.rawSensor(subjectIdx,:);
    if exist('fit405','var')
        fit405 = raw.fit405(subjectIdx,:);
    end
    dfF    = raw.dfF(subjectIdx,:);
    eventTime = round(eventTime);
    
    %behavior2 = decimate(zeros(size(photom1,1),1),dsFactor)-1.6; behavior2(eventTime) = -0.8;
    %plot(timeBehav,behavior,'b',...
    figure(figureN)
    subplot(3,1,1)
    plot(time, raw405);
    box off
    hold on
    plot(time, raw470);
    if exist('fit405','var')
        plot(time, fit405);
    end
    xline(eventTime) %only works on Matlab R2020a onwards        
    ylabel('raw fluorescence ') %                      Events per second')
    xlabel('Time (sec)')
    legend('raw405','raw470','fit405','Location','Best')
    title(['raw traces and the 405fit of ' subName{subjectIdx}])

    subplot(3,1,2)
    plot(time, 100*dfF); 
    box off
    hold on
    xline(eventTime) %only works on Matlab R2020a onwards        
    ylabel('dfF') %                      Events per second')
    xlabel('Time (sec)')
    title(['dF/F% for ' subName{subjectIdx}])

    subplot(3,1,3)
    plot(time, data1); 
    box off
    hold on
    xline(eventTime) %only works on Matlab R2020a onwards        
    ylabel('z-score') %                      Events per second')
    xlabel('Time (sec)')
    title(['raw z-score trace for ' subName{subjectIdx}])
    
    
%    yticks([-10 -5 0 5 10 15 20 25 30 35 40])
%    yticklabels({'-10','-5','0','5','10','5','10','15','20','25','30'})

end