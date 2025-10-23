function [datOUT,ts,timing] = tdtCONVERT(details, tankdir, tankname, blockname)
%  Function to extract all tdt files and data streams from their proprietary structures.
%  static refers to the static system in the Bruchas lab
%  auto refers to the automated stim protocol on the Bruchas cart

fearON = details.fearON;
stimON = details.stimON;
static = details.static;
auto = details.auto;
dual=details.dualColor;
twomice=details.twomice;

if static == 0
    storenames = {'470A'}; % name of stores to extract from TDT (usu. 4-letter code) 
    storenames2 = {'405A'};
    storenames3 = {'Epo1'};%stim
    storenames3 = {'Pe1/'};
else
    storenames = {'465A'}; % name of stores to extract from TDT (usu. 4-letter code) 
    storenames2 = {'405A'};
    storenames3 = {'Pu1/  '}; %stim
end
%LMag is the demodulated data, may also have other timestamps etc

storenames4 = {'U11_'}; % Fear 

if auto
    storenames3 = {'Pe1/'}; %{'Wi2/'};     %'Per2' or 'Wi2/'
end

if dual && ~static
    storenames = {'435A'};  % Isobestic
    storenames2 = {'490A'}; %GFP
    storenames3 = {'565B'}; %RFP
elseif dual && static
    storenames = {'405A'};  % Isobestic
    storenames2 = {'465A'}; %GFP
    storenames3 = {'560B'}; %RFP    
end

% extract
for k = 1:numel(storenames)
  if details.twomice == 1
    if contains(tankname,'Daniel') && details.twomice == 1
        storenames = {'470A'}; % name of stores to extract from TDT (usu. 4-letter code) 
        %LMag is the demodulated data, may also have other timestamps etc
    
        storenames_iso = {'405A'};
    
        storenames2     = {'470B'}; %stim
    
        storenames2_iso = {'405B'}; % Fear 
    elseif ~contains(tankname,'Daniel') && details.twomice == 1
        storenames = {'465A'}; % name of stores to extract from TDT (usu. 4-letter code) 
        %LMag is the demodulated data, may also have other timestamps etc
    
        storenames_iso = {'405A'};
    
        storenames2     = {'465C'}; %stim
    
        storenames2_iso = {'405C'}; % Fear 
    elseif details.twoGFP == 1 && details.twomice ~= 1
        storenames = {'470A'}; % name of stores to extract from TDT (usu. 4-letter code) 
        %LMag is the demodulated data, may also have other timestamps etc
    
        storenames_iso = {'405A'};
    
        storenames2     = {'470B'}; %BLA
    
        storenames2_iso = {'405B'};  
    elseif details.twoGFP == 1 && details.twomice == 1
        storenames = {'465A'}; % name of stores to extract from TDT (usu. 4-letter code) 
        %LMag is the demodulated data, may also have other timestamps etc
    
        storenames_iso = {'405A'};
    
        storenames2     = {'465C'}; %BLA
    
        storenames2_iso = {'405C'};  

    end
    
for k = 1:numel(storenames)
  storename = storenames{k};
  S{k} = tdt2mat(tankdir, tankname, blockname, storename);
  
  storenames_iso = storenames_iso{k};
  S_iso{k} = tdt2mat(tankdir, tankname, blockname, storenames_iso);

  storename2 = storenames2{k};
  S2{k} = tdt2mat(tankdir, tankname, blockname, storename2);
  
  storename2_iso = storenames2_iso{k};
  S2_iso{k} = tdt2mat(tankdir, tankname, blockname, storename2_iso);
end
    % Massage data and get time stamps
    %Fibpho2 is 465C and 405C = BoxD on left
    
    LMag = S{1}; %add more if you extracted more stores above
    % LMag2 = S{2};
    % For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
    chani1 = LMag.channels==1;
    chani2 = LMag.channels==2;
    
    LMag2 = S_iso{1}; %add more if you extracted more stores above
    % LMag2 = S{2};
    % For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
    chani21 = LMag2.channels==1;
    chani22 = LMag2.channels==2;
    
    LMag3 = S2{1}; %add more if you extracted more stores above
    % LMag3 = S3{2};
    % For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
    chani31 = LMag3.channels==1;
    chani32 = LMag3.channels==2;
    
    LMag4 = S2_iso{1}; %add more if you extracted more stores above
    % LMag3 = S3{2};
    % For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
    chani41 = LMag4.channels==1;
    chani42 = LMag4.channels==2;
    
    % Get LMag data as a vector (repeat for each channel)
    dat470D = LMag.data(chani1,:);
    dat470D = reshape(dat470D', [],1); % unwrap data from m x 256 array
    % dat405A = LMag.data(chani21,:);
    % dat405A = reshape(dat405A', [],1); % unwrap data from m x 256 array
    
    % Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
    ts = LMag.timestamps(chani1);
    t_rec_start = ts(1);
    dt = datetime( t_rec_start, 'ConvertFrom', 'posixtime' );
    
    ts = ts-ts(1); % convert from Unix time to 'seconds from block start'
    ts = bsxfun(@plus, ts(:), (0:LMag.npoints-1)*(1./LMag.sampling_rate));
    ts = reshape(ts',[],1);
    
    %%%%%%%%%%%%%%%%%% 
    dat405D = LMag2.data(chani21,:);
    dat405D = reshape(dat405D', [],1); % unwrap data from m x 256 array
    % dat405A = LMag.data(chani21,:);
    % dat405A = reshape(dat405A', [],1); % unwrap data from m x 256 array
    
    % Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
    ts2 = LMag2.timestamps(chani21);
    t_rec_start2 = ts2(1);
    
    ts2 = ts2-ts2(1); % convert from Unix time to 'seconds from block start'
    ts2 = bsxfun(@plus, ts2(:), (0:LMag2.npoints-1)*(1./LMag2.sampling_rate));
    ts2 = reshape(ts2',[],1);
    
    % 465 C and 405 C below
    dat465C = LMag3.data(chani31,:);
    dat465C = reshape(dat465C', [],1); % unwrap data from m x 256 array
    
    % Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
    ts3 = LMag3.timestamps(chani31);
    t_rec_start = ts3(1);
    dt = datetime( t_rec_start, 'ConvertFrom', 'posixtime' );
    ts3 = bsxfun(@plus, ts3(:), (0:LMag.npoints-1)*(1./LMag.sampling_rate));
    ts3 = reshape(ts3',[],1);
    
    dat405C = LMag4.data(chani41,:);
    dat405C = reshape(dat405C', [],1); % unwrap data from m x 256 array
    
    ts4 = LMag4.timestamps(chani41);
    t_rec_start4 = ts4(1);
    dt4 = datetime( t_rec_start4, 'ConvertFrom', 'posixtime' ); 
    
    ts4 = ts4-ts4(1); % convert from Unix time to 'seconds from block start'
    ts4 = bsxfun(@plus, ts4(:), (0:LMag4.npoints-1)*(1./LMag4.sampling_rate));
    ts4 = reshape(ts4',[],1);
    timingStart = seconds(dt4 - dt); 
    
    d = [size(dat470D,1),size(dat405D,1)]; [minD, minIdx] = min(d);
    c = [size(dat465C,1),size(dat405C,1)]; [minC, minIdx] = min(c);
    dat470D = dat470D(1:minD); dat405D = dat405D(1:minD);
    dat470C = dat465C(1:minC); dat405C = dat405C(1:minC);
    ts = ts(1:minD); ts2 = ts2(1:minD);
    ts3 = ts3(1:minC); ts4 = ts4(1:minC);
    dat470.dat470D = dat470D;
    dat470.dat470C = dat470C;
    dat405.dat405D = dat405D;
    dat405.dat405C = dat405C;       
 %
 elseif details.twoGFP == 1 
            storenames = {'470A'}; % FibPho1 == CA1
            storenames_iso = {'405A'};
            storenames2     = {'470B'}; %FibPho2 == BLA
            storenames2_iso = {'405B'}; % 
         if details.static == 1
                storenames = {'465A'}; % name of stores to extract from TDT (usu. 4-letter code) 
                %LMag is the demodulated data, may also have other timestamps etc
            
                storenames_iso = {'405A'};
            
                storenames2     = {'465C'}; %BLA
            
                storenames2_iso = {'405C'};  
         end
                
        for k = 1:numel(storenames)
          storename = storenames{k};
          S{k} = tdt2mat(tankdir, tankname, blockname, storename);
          
          storenames_iso = storenames_iso{k};
          S_iso{k} = tdt2mat(tankdir, tankname, blockname, storenames_iso);
        
          storename2 = storenames2{k};
          S2{k} = tdt2mat(tankdir, tankname, blockname, storename2);
          
          storename2_iso = storenames2_iso{k};
          S2_iso{k} = tdt2mat(tankdir, tankname, blockname, storename2_iso);
        end
    % Massage data and get time stamps
    
    LMag = S{1}; %add more if you extracted more stores above
    % LMag2 = S{2};
    % For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
    chani1 = LMag.channels==1;
    chani2 = LMag.channels==2;
    
    LMag2 = S_iso{1}; %add more if you extracted more stores above
    % LMag2 = S{2};
    % For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
    chani21 = LMag2.channels==1;
    chani22 = LMag2.channels==2;
    
    LMag3 = S2{1}; %add more if you extracted more stores above
    % LMag3 = S3{2};
    % For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
    chani31 = LMag3.channels==1;
    chani32 = LMag3.channels==2;
    
    LMag4 = S2_iso{1}; %add more if you extracted more stores above
    % LMag3 = S3{2};
    % For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
    chani41 = LMag4.channels==1;
    chani42 = LMag4.channels==2;
    
    if fearON
        storename4 = storenames4{k}
        S5{k} = tdt2mat(tankdir, tankname, blockname, storename4);
        LMag5 = S5{1}; %add more if you extracted more stores above
        % LMag3 = S3{2};
        % For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
        chani51 = LMag5.channels==1;
        chani52 = LMag5.channels==2;
    end


    % Get LMag data as a vector (repeat for each channel)
    dat470A = LMag.data(chani1,:);
    dat470A = reshape(dat470A', [],1); % unwrap data from m x 256 array
    % dat405A = LMag.data(chani21,:);
    % dat405A = reshape(dat405A', [],1); % unwrap data from m x 256 array
    
    % Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
    ts = LMag.timestamps(chani1);
    t_rec_start = ts(1);
    dt = datetime( t_rec_start, 'ConvertFrom', 'posixtime' );
    
    ts = ts-ts(1); % convert from Unix time to 'seconds from block start'
    ts = bsxfun(@plus, ts(:), (0:LMag.npoints-1)*(1./LMag.sampling_rate));
    ts = reshape(ts',[],1);
    
    %%%%%%%%%%%%%%%%%%
     
    dat405A = LMag2.data(chani21,:);
    dat405A = reshape(dat405A', [],1); % unwrap data from m x 256 array
    % dat405A = LMag.data(chani21,:);
    % dat405A = reshape(dat405A', [],1); % unwrap data from m x 256 array
    
    % Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
    ts2 = LMag2.timestamps(chani21);
    t_rec_start2 = ts2(1);
    
    ts2 = ts2-ts2(1); % convert from Unix time to 'seconds from block start'
    ts2 = bsxfun(@plus, ts2(:), (0:LMag2.npoints-1)*(1./LMag2.sampling_rate));
    ts2 = reshape(ts2',[],1);
    
    % 465 C and 405 C below
    dat470B = LMag3.data(chani31,:);
    dat470B = reshape(dat470B', [],1); % unwrap data from m x 256 array
    
    % Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
    ts3 = LMag3.timestamps(chani31);
    t_rec_start = ts3(1);
    dt = datetime( t_rec_start, 'ConvertFrom', 'posixtime' );
    ts3 = bsxfun(@plus, ts3(:), (0:LMag.npoints-1)*(1./LMag.sampling_rate));
    ts3 = reshape(ts3',[],1);
    
    dat405B = LMag4.data(chani41,:);
    dat405B = reshape(dat405B', [],1); % unwrap data from m x 256 array
    
    ts4 = LMag4.timestamps(chani41);
    t_rec_start4 = ts4(1);
    dt4 = datetime( t_rec_start4, 'ConvertFrom', 'posixtime' ); 
    
    ts4 = ts4-ts4(1); % convert from Unix time to 'seconds from block start'
    ts4 = bsxfun(@plus, ts4(:), (0:LMag4.npoints-1)*(1./LMag4.sampling_rate));
    ts4 = reshape(ts4',[],1);
    timingStart = seconds(dt4 - dt); 
    
    a = [size(dat470A,1),size(dat405A,1)]; [minA, minIdx] = min(a);
    b = [size(dat470B,1),size(dat405B,1)]; [minB, minIdx] = min(b);
    dat470A = dat470A(1:minA); dat405A = dat405A(1:minA);
    dat470B = dat470B(1:minB); dat405B = dat405B(1:minB);
    ts = ts(1:minA); ts2 = ts2(1:minA);
    ts3 = ts3(1:minB); ts4 = ts4(1:minB);
    dat470.dat470A = dat470A;
    dat470.dat470B = dat470B;
    dat405.dat405A = dat405A;
    dat405.dat405B = dat405B;       


    
  else
      storename = storenames{k};
          S{k} = tdt2mat(tankdir, tankname, blockname, storename);
          
          storename2 = storenames2{k};
          S2{k} = tdt2mat(tankdir, tankname, blockname, storename2);
        if stimON || dual
          storename3 = storenames3{k};
          S3{k} = tdt2mat(tankdir, tankname, blockname, storename3);
        end
        if fearON
          storename4 = storenames4{k}
          S5{k} = tdt2mat(tankdir, tankname, blockname, storename4);
        end

        % Massage data and get time stamps
        
        LMag = S{1}; %add more if you extracted more stores above
        % LMag2 = S{2};
        % For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
        chani1 = LMag.channels==1;
        chani2 = LMag.channels==2;
        
        LMag2 = S2{1}; %add more if you extracted more stores above
        % LMag2 = S{2};
        % For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
        chani21 = LMag2.channels==1;
        chani22 = LMag2.channels==2;
        % chani21 = LMag2.channels==1;
        % chani22 = LMag2.channels==2;
        
        if stimON || dual
        LMag3 = S3{1}; %add more if you extracted more stores above
        % LMag3 = S3{2};
        % For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
        chani31 = LMag3.channels==1;
        chani32 = LMag3.channels==2;
        end
        
        if fearON
        LMag5 = S5{1}; %add more if you extracted more stores above
        % LMag3 = S3{2};
        % For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
        chani51 = LMag5.channels==1;
        chani52 = LMag5.channels==2;
        end
        
        %{
        % Get LMag data as a vector (repeat for each channel)
        if contains(LMag.storename,'405')
            dat405 = LMag.data(chani1,:);
            dat405 = reshape(dat405', [],1); % unwrap data from m x 256 array
            
            dat470 = LMag2.data(chani21,:);
            dat470 = reshape(dat470', [],1); % unwrap data from m x 256 array
        else
            warning('the 405 and 470 were switched')
            dat470 = LMag.data(chani1,:);
            dat470 = reshape(dat470', [],1); % unwrap data from m x 256 array
            
            dat405 = LMag2.data(chani21,:);
            dat405 = reshape(dat405', [],1); % unwrap data from m x 256 array
        end
        %}

        % dat405 = LMag.data(chani21,:);
        % dat405 = reshape(dat405', [],1); % unwrap data from m x 256 array
        
        % Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
        ts = LMag.timestamps(chani1);
        t_rec_start = ts(1);
        
        ts = ts-ts(1); % convert from Unix time to 'seconds from block start'
        ts = bsxfun(@plus, ts(:), (0:LMag.npoints-1)*(1./LMag.sampling_rate));
        ts = reshape(ts',[],1);
        
        % dat405 = LMag.data(chani21,:);
        % dat405 = reshape(dat405', [],1); % unwrap data from m x 256 array
        
        % Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
        ts2 = LMag2.timestamps(chani21);
        t_rec_start2 = ts2(1);
        dt = datetime( t_rec_start2, 'ConvertFrom', 'posixtime' );
        
        ts2 = ts2-ts2(1); % convert from Unix time to 'seconds from block start'
        ts2 = bsxfun(@plus, ts2(:), (0:LMag2.npoints-1)*(1./LMag2.sampling_rate));
        ts2 = reshape(ts2',[],1);

        if dual
            dat565 = LMag3.data(chani31,:);
            dat565 = reshape(dat565', [],1); % unwrap data from m x 256 array            
        end



            % Get LMag data as a vector (repeat for each channel)
        dat470A = LMag.data(chani1,:);
        dat470A = reshape(dat470A', [],1); % unwrap data from m x 256 array
        % dat405A = LMag.data(chani21,:);
        % dat405A = reshape(dat405A', [],1); % unwrap data from m x 256 array
        
        % Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
        ts = LMag.timestamps(chani1);
        t_rec_start = ts(1);
        dt = datetime( t_rec_start, 'ConvertFrom', 'posixtime' );
        
        ts = ts-ts(1); % convert from Unix time to 'seconds from block start'
        ts = bsxfun(@plus, ts(:), (0:LMag.npoints-1)*(1./LMag.sampling_rate));
        ts = reshape(ts',[],1);
        
        %%%%%%%%%%%%%%%%%%
         
        dat405A = LMag2.data(chani21,:);
        dat405A = reshape(dat405A', [],1); % unwrap data from m x 256 array
        % dat405A = LMag.data(chani21,:);
        % dat405A = reshape(dat405A', [],1); % unwrap data from m x 256 array
        
        % Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
        ts2 = LMag2.timestamps(chani21);
        t_rec_start2 = ts2(1);
        
        ts2 = ts2-ts2(1); % convert from Unix time to 'seconds from block start'
        ts2 = bsxfun(@plus, ts2(:), (0:LMag2.npoints-1)*(1./LMag2.sampling_rate));
        ts2 = reshape(ts2',[],1);
     
        a = [size(dat470A,1),size(dat405A,1)]; [minA, minIdx] = min(a);
        dat470A = dat470A(1:minA); dat405A = dat405A(1:minA);
        ts = ts(1:minA); ts2 = ts2(1:minA);
   
        dat470.dat470A = dat470A;

        dat405.dat405A = dat405A;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For stim and Fear 

if dual

    storenames = {'435A'};  % Isobestic
    storenames2 = {'490A'}; %GFP
    storenames3 = {'565B'}; %RFP
    datOUT.dat435 = dat405;
    datOUT.dat490 = dat470;
    datOUT.dat565 = dat565;
elseif details.twoGFP
    datOUT.dat470A = dat470A;
    datOUT.dat470B = dat470B;
    datOUT.dat405A = dat405A;
    datOUT.dat405B = dat405B;
else
    datOUT.dat470 = dat470A;
    datOUT.dat405 = dat405A;
end

if stimON
% Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
ts3 = LMag3.timestamps(chani31);

%multiple stims
if any(diff(ts3)>1)
    dTs3 = find(diff(ts3)>1);
    stimLocs = 1;
    for i = 1:length(dTs3)
        stimLocs = [stimLocs dTs3(i)+1];
    end
end
if exist('stimLocs')
    t_rec_start3 = ts3(1);
    dt3 = datetime( t_rec_start3, 'ConvertFrom', 'posixtime' );
    ts3St = ts3(1); % start of ts3
    ts3 = ts3-ts3(1); % convert from Unix time to 'seconds from block start'
    timing = seconds(dt3 - dt);
    for i = 2:length(stimLocs)
        timing = [timing timing(1)+ts3(stimLocs(i))];
    end
else
    t_rec_start3 = ts3(1);
    dt3 = datetime( t_rec_start3, 'ConvertFrom', 'posixtime' );
    ts3St = ts3(1); % start of ts3
    ts3 = ts3-ts3(1); % convert from Unix time to 'seconds from block start'
    timing = seconds(dt3 - dt);
end
end
end


if fearON 
dat5 = LMag5.data(chani51,:);
dat5 = reshape(dat5', [],1); % unwrap data from m x 256 array

% Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
ts5 = LMag5.timestamps(chani51);
t_rec_start5 = ts5(1:length(ts5)); t_rec_start5=t_rec_start5(logical(dat5));
dt5 = datetime( t_rec_start5, 'ConvertFrom', 'posixtime' ); 

ts5 = ts5-ts5(1); % convert from Unix time to 'seconds from block start'
ts5 = bsxfun(@plus, ts5(:), (0:LMag5.npoints-1)*(1./LMag5.sampling_rate));
ts5 = reshape(ts5',[],1);
timing = seconds(dt5 - dt); 

end

if fearON && stimON
    warning('timing was overwritten, you did fearTTLs and stim!')
    warning('write more code')
end





if ~exist('timing')
    timing = [];
end
end