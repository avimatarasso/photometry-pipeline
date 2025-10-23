function [dataFilt,fitIso,f2,fitcurve,dataFiltExp] = fitAndSub(datSensor, datIso, Dts,details, SensorStr, IsoStr,subject)
%  File to convert all TDT recordings to df/F
% Apply a linear least squared fit of the isobestic to the sensor, and
% subtract the sensor signal from the isobestic. The function also plots 
% these time series

BLlength = details.BLlength;
Fs   = details.Fs;

% preprocessing, fit to linear least squares/exponential, subtract isobestic
if details.subON
    lenSensor = length(datSensor);
    lenIso = length(datIso);
    if lenSensor>lenIso %470 greater than 405
        datSensor = datSensor(1:lenIso);
        Dts = Dts(1:lenIso);
        fprintf(['Shortened sensor to iso length\n'])
    elseif lenIso>lenSensor %405 greater than 470
        datIso = datIso(1:lenSensor);
        Dts = Dts(1:lenSensor);
        fprintf(['Shortened iso to sensor length\n'])
    elseif Dts>lenSensor
        Dts = Dts(1:lenSensor);
    end    

    fitIso = LLS(datSensor,datIso);
    if details.subFit
        dataFilt = (datSensor - fitIso)./fitIso; %470 - fit(405) - F0 then normalized to 405fit
        dataFilt = dataFilt - mean(dataFilt(1:BLlength*Fs)); % Center your dF/F? CUSTOMIZE
    else
        datSensor = datSensor - datIso;  
        dataFilt = datSensor/mean(datSensor(1:BLlength*Fs));
        dataFilt = dataFilt - mean(dataFilt(1:BLlength*Fs)); % Center your dF/F? CUSTOMIZE
    end
    
    if details.check ==1
        % plot raw signals, linear least squares fit, and the df/f%
        if details.allfigures 
            figure
        else
            figure(1), clf
        end
        
        if Dts>lenSensor
            Dts = Dts(1:length(datSensor));
        elseif Dts<lenSensor
            Dts = 1/Fs:1/Fs:length(datSensor)/Fs;
            Dts = Dts';

        end  
        subplot(3,1,1)
        plot(Dts,datSensor), hold on
        plot(Dts,datIso)
        legend(SensorStr,IsoStr);%'LLSfit405', 'df/f + mean(datSensor)')
        ylim([220 280])
        title('raw traces')
        ylabel('mV')
        xlabel('time (s)')
        hold off

        subplot(3,1,2)
        plot(Dts,datSensor), hold on
        plot(Dts,fitIso)
        legend(SensorStr,['llsfit(' IsoStr ')']);%'LLSfit405', 'df/f + mean(datSensor)')
        %ylim([mean(fitIso)-30 mean(fitIso)+30])
        title([subject ' linear least squares fit - just LLS subtraction'])
        ylabel('mV')
        xlabel('time (s)')
        hold off
        
        subplot(3,1,3)
        plot(Dts,dataFilt*100)
        ylabel('dF/F %')
        xlabel('time (s)')
        ylim([-5 10])
        title(['dF/F% of ' SensorStr])
        hold off
    end
    
    % Do you want to subtract an exponential decay or LLS? 
    % Make sure the one you want is ending up in the dfF variable
    f2 = fit(Dts,datSensor,'exp2');
    fitcurve= f2(Dts);
    dataFiltExp =  datSensor - fitcurve;
    
    dataFiltExp = (datSensor - datIso)./datIso; %470 - fit(405) - F0 then normalized to 405fit
    dataFiltExp = dataFiltExp - mean(dataFiltExp(1:BLlength*Fs)); % Center your dF/F? CUSTOMIZE

    if details.check ==1
        % plot raw signals, linear least squares fit, and the df/f%
        if details.allfigures 
            figure
        else
            figure(1), clf
        end
        subplot(3,1,1)
        plot(Dts,datSensor), hold on
        plot(Dts,datIso)
        legend(SensorStr,IsoStr);%'LLSfit405', 'df/f + mean(datSensor)')
        ylim([220 280])
        title('raw traces')
        ylabel('mV')
        xlabel('time (s)')
        hold off

        subplot(3,1,2)
        plot(Dts,datSensor), hold on
        plot(Dts,fitIso)
        legend(SensorStr,[SensorStr ' expfit']);%'LLSfit405', 'df/f + mean(datSensor)')
        ylim([mean(fitIso)-30 mean(fitIso)+30])
        title([subject SensorStr ' exponential fit'])
        ylabel('mV')
        xlabel('time (s)')
        hold off
        
        subplot(3,1,3)
        plot(Dts,dataFiltExp*100)
        ylabel('dF/F %')
        xlabel('time (s)')
        ylim([-5 10])
        title(['dF/F% of (' SensorStr ' - ' IsoStr ')/' IsoStr])
        hold off
    end
    
    
    
else %if not subtracting the linear least squares, just center the data
    %fitIso = LLS(datSensor,datIso);
    %datSensor = datSensor - datIso;
    dataFilt = (datSensor - mean(dataFilt(1:BLlength*Fs)))./mean(datSensor(1:BLlength*Fs));% - fitIso)./fitIso; %470 - fit(405) - F0 then normalized to 405fit
    fitIso = nan; f2 = nan; fitcurve=nan; dataFiltExp =nan;
    if details.check ==1
        % plot raw signals, linear least squares fit, and the df/f%
        if details.allfigures 
            figure
        else
            figure(1), clf
        end
        
        plot(Dts,dataFilt*100)
        ylabel('dF/F %')
        xlabel('time (s)')
        ylim([-5 10])
        title(['dF/F% of ' SensorStr ', without ' IsoStr ' subtraction'])
        hold off
    end
    
    
end
end