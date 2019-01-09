%% Cage 8 MSCEEG2

    %% Identify optimal 3 min baseline

    tic
    timerVal = tic;

    delete(gcp('nocreate')); % start parallel pool
    p = parpool('local');
    mr = mapreducer(p);
    % mapreducer(0);

    %% Setup datastore for cage 8 data
    
    disp('cage 8 baseline')
    ds = tabularTextDatastore('10262015_A9-16_data.txt', 'TreatAsMissing', 'AMPSAT'); % txt file to analyze
    ds.SelectedVariableNames = {'Time', 'x8', 'x6', 'x31', 'x11'}; % time and MSCEEG2 cage 8

    %% reads for initial baseline

    ds.ReadSize = 32160;
    EEG1 = table2array(read(ds));
    EEG2 = table2array(read(ds));
    EEG3 = table2array(read(ds));
    EEG4 = table2array(read(ds));
    EEG5 = table2array(read(ds));
    EEG6 = table2array(read(ds));

    baseline_EEG = [EEG1;EEG2;EEG3;EEG4;EEG5;EEG6];

    %% set up variables

    best_baseline_EEG = zeros((length(baseline_EEG)),5);
    best_baseline_spread = inf;
    best_baseline_time = inf;

%% loop through potential baselines

while hasdata(ds)
    
    baseline_spread = mean(median(baseline_EEG(:,2:5))+std(baseline_EEG(:,2:5)));
    
    if baseline_spread < best_baseline_spread
        best_baseline_EEG = baseline_EEG(:,2:5);
        best_baseline_spread = baseline_spread;
        best_baseline_time = baseline_EEG(1,1);
    end
    
    EEG_step = table2array(read(ds)); % read new step in EEG
    step = length(EEG_step);
    
    baseline_EEG = [baseline_EEG(step:end,:); EEG_step]; % new baseline
    time_end = baseline_EEG(end,1);
    
end

best_baseline_EEG_8 = best_baseline_EEG;
best_baseline_spread_8 = best_baseline_spread;
best_baseline_time_8 = best_baseline_time;
time_end_8 = time_end;

clear best_baseline_EEG best_baseline_spread best_baseline_time time_end

%% Bergstrom code rearranged and edited
    
    %For UW Madison Data: 
    f = 1024; %sampling rate is 1024 Hz
    CH = 4; %update as needed
    win0 = 256; %250 ms = 1024/4
    DC = 4; %need to also try 4?
    
    
%% Bergstrom baseline code

EEG = best_baseline_EEG_8; 
    % specify files based on readout shown. 
    % to select sequential files, format x:y 
    % -or- put in a array like [a b c e]

    % adjust window size so it can be evenly divided by the downsampling
    % factor identified by the decomposition level
    dc = 2^DC;
    win = floor(win0/dc)*dc;
    w = win/dc;
    
%ww=1; % is used to set up the Decomp Variable - WTF do I need to set this for?


% Set up variables for finding the threshold, h is one number, so length=1.
% BaseStatMedian = zeros(1,length(h));
% BaseStatStDev = zeros(1,length(h));

clear win0


 
% Find median, standard deviation of baseline files. These values will be 
% used to set the threshold for event identification. 
% 
% 
% 
% Perform this loop on baseline file(s) only 
% 
% h = baseline file id (from input above)

%%
for k = 1   
    
%%
        %load (files(k).name); % load the k/h-th (baseline) file in files 
        %EEG = clean; % clean is the EEG variable, all channels.
  %%
  
        for col = 1:CH    
                    % perform this loop on each channel individually
   %%     
            [C,L]  = wavedec(EEG(:,col),DC,'db4');
                    % perform wavelet decomposition on the designated channel
                    % to the indicated decomposition level (DC) using the 
                    % db4 wavelet
                    % for DC = 4 on 400 hz data, pull out XXXXX frequency
                    % for DC = 5 on 1024 hz data, pull out yyyyy frequency
                    
                    % From help wavedec
                    % The output vector, C, contains the wavelet decomposition. 
                    % L contains the number of coefficients by level.
                    % C and L are organized as:
                    % C      = [app. coef.(N)|det. coef.(N)|... |det. coef.(1)]
                    % L(1)   = length of app. coef.(N)
                    % L(i)   = length of det. coef.(N-i+2) for i = 2,...,N+1
                    % L(N+2) = length(X).
               
               
            %L2 = [0; L]; % column of L with 0 added to the beginning
         %   Decomp = C((sum(L2(1:ww-1))+1):L(ww)); 
         Decomp = C(1:L(1));
                [rows,~] = size (Decomp); % '~' voids the 1 column. 
                    % Set Decomp equal to the approximation at the DC level
                    % The first value of L gives the number of coefficients 
                    % from the beginning of C that contains the 4th level
                    % approximation from the transform. This contains the 
                    % "most information" while reducing noise.


  %%                  
                    
                    % find distance between datapoints n and n-1
        Dist = zeros(rows,1); 
            % creates a column vectors full of zeros that is the size of 
            % the decomposed/approximation of the EEG (rows).
            for m = 2:(rows-w) % size of Decomp (rows) - window (win = 96)
                %w = win/16 or win/DC^4
                
                Dist((m-1)) = abs(Decomp(m)-Decomp(m-1));
                    % modifies Dist and inputs the individual line lengths. 
                    % Starts from the 1st row in Dist and gives the
                    % absolute distance between the 2nd point and 1st point
                    % of the approximation/decomposed EEG.
                    % Ends at rows-6-1 and gives the abs diff between
                    % rows-6 and rows-6-1
            end
     %%                              
                    % sum distance between datapoints across window
           % LLength = zeros(length(Dist),col); 
                % creates LLength full of a column of zeros that has a
                % length of Dist which is the size of the decomp EEG
            LLStart = sum(Dist(1:w)); 
                % sums up line lengths of first window (6 numbers long,
                % represents ~100 datapoints from window)
                % which is Dist 1 to Dist 6,
                % which are line lengths between points 1 to 97 of the 
                % decomp EEG
            LLength(1) = LLStart;
                % modifies the first entry of LLength to the sum of line
                % lengths of first window
  %%              
        for y = 2:rows-w 
            % this is the new LLength entry row, aka the 2nd window
                % which is sum of Dist 2 to Dist 97
                % which are line lengths between points 2 to 98 of the
                % decomp EEG
            % same as m, the number of times we calculated distance
                % runs from from 2 to rows-window(96) 
                % PROBLEM: the last window, from rows-96 to rows has 0 dist
                % because it wasn't calculated in Dist above.
            yy = y+w-1; 
                % starts at 2+96-1=97 to rows-96+96-1=rows-1
                    % specifies the new Dist entry to add to this 
                    % new window.
            LL_a = LLength(y-1)-Dist(y-1); 
                % LLength(2-1) - Dist(2-1),
                % LLength(rows-96-1) - Dist(rows-96-1)   
                    % takes the LL sum of previous window and subtracts
                    % the first Dist entry of that previous window
            LLength(y) = LL_a + Dist(yy);
                % Adds the next Dist entry to the new window 
                % Modifies LLength with this new sum of dist in new window
        end
        %%
        
        j = 1;
            % j = baseline file number in 'files' 
            % - first baseline file number in 'files' 
            % + 1, answer is 1 if there is only one baseline file 
            % used if there are multiple baseline files.
                    
                    % find median and standard deviation of calculated
                    % line lengths
       
            BaseStatMedian(j,col) = median(LLength);
            BaseStatStDev(j,col) = std(LLength);
            
            clear Decomp rows Dist LLength LLStart LL_a y yy C L
        end
%%
    clear EEG 
 
%%
    
end

%if length(h)>1 
        % the stats are averaged across baseline files and logged into a
        % BaseStats variable.
%     BaseStats(1,:) = BaseStatMedian;
%     BaseStats(2,:) = BaseStatStDev;
    
     BaseStats = [BaseStatMedian; BaseStatStDev];
%else
 %   BaseStats(1,:) = BaseStatMedian;
  %  BaseStats(2,:) = BaseStatStDev;
%end

%clear BaseStatMedian BaseStatStDev j


%%

SF = 4.0;       % SF can be adjusted to increase or decrease the threshold 
                % used to determine if a line length is different from the 
                % median baseline line length. SF is a multiplier and
                % identifies how many standard deviations from the median a
                % line length value must be in order for it to be pulled
                % out in the analysis.

%% Setup datastore for detections

ds = tabularTextDatastore('10262015_A9-16_data.txt', 'TreatAsMissing', 'AMPSAT');
ds.SelectedVariableNames = {'x8', 'x6', 'x31', 'x11', 'Time'}; % time and MSCEEG2 cage 8
ds.ReadSize = 100000;

total_seizures = []; % setup variables for detections
total_abnormal = [];
total_spikes = [];

%% Detection pass 1

disp('seizure detection: pass 1')

while hasdata(ds)
    
    EEG_1 = table2array(read(ds));
    if hasdata(ds) EEG_2 = table2array(read(ds)); else EEG_2 = []; end;
    EEG = [EEG_1; EEG_2];
    progress = 100*(EEG(1,5)/time_end_8);
    elapsedTime = toc(timerVal);
    fprintf('Progress %0.1f %% ', progress)
    fprintf('in %0.1f min \r', elapsedTime/60)
    
for k = 1
%%
  
  
%     load (files(k).name);
%         input_name = files(k).name;            
%         [~, name, ~] = fileparts (input_name); 
        %EEG = clean;  
            fr=length(EEG)/f;
            t=EEG(:,5);             
            %t=t(:); 
    %%
    for col = 1:CH  % run all channels
                    
        [C,L]  = wavedec(EEG(:,col),DC,'db4');
          % perform wavelet decomposition to the indicated 
                    % decomposition level (DC) using the db4 wavelet
            L2 = [0; L];
             Decomp = C(1:L(1));
        %Decomp = C((sum(L2(1:ww-1))+1):L(ww));
            [rows, ~] = size (Decomp);
            % Set Decomp equal to the approximation at the DC level
       %%  
            % find distance between datapoints n and n-1
        Dist = zeros(rows,1);
            for r = 2:rows-w
                Dist(r-1) = abs(Decomp(r)-Decomp(r-1));
            end
            %%
            % sum distance between datapoints across window
            LLength = zeros(rows,1);
            LLStart = sum(Dist(1:w));
            LLength(1) = LLStart;
                for r = 2:rows-w
                    LL_a = LLength(r-1)-Dist(r-1);
                    LLength(r) = LL_a+Dist(r+w);
                end
            LLAve = LLength-BaseStats(1,col);
        
            clear rows LLAveSD
            [rows, ~] = size(Dist);
      
            LLAveSD = zeros(rows,1);
            LLAveSD = LLAve/(SF*BaseStats(2,col)); 
         
        % create binary vector of thresholded data
            [rows, columns] = size(LLAveSD);
                LLTh = zeros(rows,1);
                for r = 1:rows
                    if LLAveSD(r) >= 1
                        LLTh(r) = 1;
                    else, LLTh(r) = 0;
                    end
                end
                
        % erode (open) data to remove noise - half-window size = 1.
            LLTh2 = zeros(size(LLTh));
            nhood = 1;
            GG = length(LLTh);
                for r = 1
                    if sum(LLTh(1:(r+nhood))) == length(LLTh(1:(r+nhood)))
                        LLTh2(r) = 1;    
                    else, LLTh2(r) = 0;
                    end
                end
                for r = 2:GG-1
                    if sum(LLTh(r-nhood:r+nhood)) == nhood*2+1
                        LLTh2(r)=1;      
                    else, LLTh2(r) = 0;
                    end
                end
                for r = GG-1:GG
                    if sum(LLTh(r-nhood:GG)) == length(LLTh((r-nhood:GG)))
                        LLTh2(r) = 1;         
                    else, LLTh2(r) = 0;
                    end
                end
            
        % Dilate data to restore eroded points - half-window size = 1
            LLTh3 = zeros(size(LLTh));
                for r = 1
                    if sum(LLTh2(1:(r+nhood))) >= 1
                        LLTh3(r) = 1;      
                    else, LLTh3(r) = 0;
                    end
                end
                for r = 2:GG-1
                    if sum(LLTh2(r-nhood:r+nhood)) >=1
                        LLTh3(r)=1;          
                    else, LLTh3(r) = 0;
                    end
                end
                for r = GG-1:GG
                    if sum(LLTh2(r-nhood:GG)) >= 1
                        LLTh3(r) = 1; 
                    else, LLTh3(r) = 0;
                    end
                end
        
        % Dyadic upsampling at odd indices must be performed once for each
        % decomposition level.
            dse0 = LLTh3;
            dse1 = dyadup(dse0,1); %dyadic upsampling, odd indices (1)
            dse2 = dyadup(dse1,1); %dyadic upsampling, odd indices (2)
            dse3 = dyadup(dse2,1); %dyadic upsampling, odd indices (3)
            dse4 = dyadup(dse3,1); %dyadic upsampling, odd indices (4)
            dse5 = dyadup(dse4,1); 

           % if DC = 4,
           %     dse_U = dse4;
           % elseif DC = 5
           %     dse_U = dse5;
           % end
            
        % close events - half-window = half * decomp level squared to close
        % the effect of dyadic upsampling
            GG = length(dse4);
            nhood = 0.5*dc; 
            dse4_c = zeros(GG,2);
                for r = 1:nhood
                    if sum(dse4(1:(r+nhood))) ~= 0
                        dse4_c(r) = 1;
                        
                    else, dse4_c(r) = 0;
                    end
                end
                for r = nhood+1:GG-nhood
                    if sum(dse4(r-nhood:r+nhood)) ~=0
                        dse4_c(r)=1;
                    else, dse4_c(r) = 0;
                    end
                end
                for r = GG-nhood+1:GG
                    if sum(dse4(r-nhood:GG)) ~=0
                        dse4_c(r)=1;
                        
                    else, dse4_c(r) = 0;
                    end
                end

        % Erode the ends of events
                dse4_f = zeros(size(dse4_c));
                for r = 1:nhood
                    if sum(dse4_c(1:(r+nhood))) == ...
                            length(dse4_c(1:(r+nhood)))
                        dse4_f(r) = 1;
                    else ,dse4_f(r) = 0;
                    end
                end
                for r = nhood+1:GG-nhood
                    if sum(dse4_c(r-nhood:r+nhood)) == nhood*2+1
                        dse4_f(r)=1;
                    else ,dse4_f(r) = 0;
                    end
                end
                for r = GG-nhood+1:GG
                    if sum(dse4_c(r-nhood:GG)) == ...
                            length(dse4_c((r-nhood):GG))
                        dse4_f(r) = 1;
                    else ,dse4_f(r) = 0;
                    end
                end
        
            LL_above_th_row = find(dse4_f==1); % edits to convert row to time 2/2 datastore
            LL_above_th = EEG(LL_above_th_row,5);    
                
                % pick out line lengths greater than 
                % BaseMedian + SF*StDev. return as 
                % row index (time) of event start
                       
            clear LLAve LLAveSD LLTh LLTh2 LLTh3 nhood GG r dse0 dse1 ...
                dse2 dsde3 dse4 dse4_c 
 %%               
    if isempty(LL_above_th)~=1
        
        %% Event classification    
        
        % RJK event times
        
        % organize all threshold crossings with start, end, length
        all_above_th = zeros((length(LL_above_th)-1), 3);           
        all_above_th(:,1) = LL_above_th(1:(length(LL_above_th)-1),1);
        all_above_th(:,2) = LL_above_th(2:end,1);
        all_above_th(:,3) = all_above_th(:,2) - all_above_th(:,1);
        
        % identify and organize events as those with > 750ms intervals
        % between
        if size(all_above_th,1) == 0, events = [0 0 0];
        elseif size(all_above_th,1) == 1, events = all_above_th;
        else events = zeros(length(find(all_above_th(:,3) > .75)),3); ...   
            events(:,2) = all_above_th(find(all_above_th(:,3) > 0.75),1); ...
            events(:,1) = all_above_th((find(all_above_th(:,3) > 0.75)),2); ...
            events(:,1) = circshift(events(:,1),1); ...
            events(1,1) = all_above_th(1,1);
        end
        
        % if size(events,1) >1 events(2:end,1) = events(1:(size(events,1)-1),1); ...
        %    elseif size(events,1) ~= 0 events(1,1) = all_above_th(1,1);, else events = [0 0 0];,end;
        %if length(events) ~= 0 events(:,3) = events(:,2) - events (:,1);, else events = [];,end;
        
            % identify hits that are greater than 1 from neighboring event
%                evnt_start = zeros(length(LL_above_th),1);  
%                     for r = 2:length(LL_above_th)
%                         evnt_start(r,1) = LL_above_th(r) - LL_above_th(r-1);
%                     end
%                 evnt_indx = find(evnt_start > 0.001); % edited for actual time vs rows
%        
%         % Find event lengths
%                 evnt_lngth = zeros(length(evnt_indx),1);
%                     for r = 1:length(evnt_indx)-1
%                         evnt_lngth(r,1) = LL_above_th(evnt_indx(r+1)-1) ...
%                             - LL_above_th(evnt_indx(r));
%                     end
%                 evnt_lngth(length(evnt_indx),1) = LL_above_th(end,1)...
%                     - LL_above_th(evnt_indx(end,1));
% 
%                 [row1, col1] = size(evnt_lngth);
% 
%         % Find event ends
%                 evnt_end = zeros(size(evnt_indx));
%                     for r = 1:length(evnt_indx-1)
%                         evnt_end(r,1) = LL_above_th(evnt_indx(r))+ evnt_lngth(r);
%                     end
%                 ev_end = [LL_above_th(evnt_indx) evnt_end];
%               %  ev_end(1,3) = 400;
%                 ev_end(1,3) = f;
%                     for r = 2:length(evnt_indx)
%                             ev_end(r,3) = ev_end((r),1)-ev_end((r-1),2);
%                     end
%     
%         % Identify events less than 750 ms from neighbor.
%         % Events closer than 750 ms to neighbor will be combined with neighbor.
%                 ev_end(:,4) = 0;
%                     for r=1:length(evnt_indx)
%                         if ev_end(r,3) <= .75*f
%                            ev_end(r,4) = 0;
%                         else ,ev_end(r,4) = 1;
%                         end
%                     end
%                 n_to_combine = find(ev_end(:,4) == 1); 
%                     % index of events closer than 250 ms
%     
%         % recreate event list after combining neighbors 
%                 events = zeros(length(n_to_combine),2);
%                     for r = 1:length(n_to_combine)-1
%                         events(r,1) = ev_end(n_to_combine(r),1);
%                         events(r,2) = ev_end((n_to_combine(r+1)-1),2);
%                     end
%                 events(length(n_to_combine),1) = ...
%                     ev_end(n_to_combine(end),1);
%                 events(length(n_to_combine),2) = ev_end(end,2);
    
        % Sort events based on event lengths
            % Short = events < 3 sec, includes spikes and abnormal events; 
            % Seizures = events > 3 sec in length 
                events(:,3) = events(:,2)-events(:,1);
                [rows, columns] = size(events);
                    for r = 1:rows
                        %if events(r,3) < 5*f
                         if events(r,3) < 3 % adjusted from 1 (sec?) to 3
                            events(r,4) = 1;
                        else ,events(r,4) = 2;
                        end
                    end
                ev1 = find(events(:,4)==1);
                ev2 = find(events(:,4)==2);

                    
                    if isempty(ev2) ~= 1
                        seizures = events(ev2,1:2); % removed division by f
                    else ,seizures = [];
                    end

        % Sort Short events based on event amplitude
            % Spikes = events > 250 uV in amplitude
            % Abnormal = events < 250 uV in amplitude
                if isempty(ev1) ~= 1 
                    short = events(ev1,1:2);
                    [rows, ~] = size(short);
                    SpMin = zeros(rows,2);
                    
                   EEGev = zeros(length(dse4_f),2);
                   EEGev(1:length(EEG),2) = EEG(:,5); % added time column to EEGev
                   EEGev(1:length(EEG),1) = EEG(:,col);
                    for r = 1:rows
                        sp = EEGev(find(short(r,1)==EEGev(:,2)):find(short(r,2)==EEGev(:,2)),1);
                        % sp = sp';
                        if sp ~=0 SpMinMax(r,:) = [min(sp) max(sp)];, else SpMinMax = zeros(r,2); end;
                    end
        
                SpAmp = SpMinMax(:,2)-SpMinMax(:,1);
                    Abn0 = find(SpAmp <.250);
                        abnormal = short(Abn0,1:2);
                    Sp0 = find(SpAmp>=.250);
                        spikes = short(Sp0,1:2);    
                else
                        spikes = [];
                        abnormal = [];
                end

                spikecount = length(spikes);
                abnormalcount = length(abnormal);
                seizurecount = length(seizures);
                    counts = [spikecount seizurecount abnormalcount ];
          
     % find the actual abnormal time totals
    
                if isempty(spikes) ~=1
                    TotSPTime_byAlg = sum(spikes(:,2)-spikes(:,1));
                else ,TotSPTime_byAlg = 0;
                end
                if isempty(abnormal)~=1
                    TotABTime_byAlg = sum(abnormal(:,2)-abnormal(:,1));
                else ,TotABTime_byAlg = 0;
                end
                if isempty(seizures) ~= 1
                    TotSZTime_byAlg = sum(seizures(:,2)-seizures(:,1));
                else ,TotSZTime_byAlg = 0;
                end
                    RealTotTime = [TotSPTime_byAlg  ...
                        TotSZTime_byAlg TotABTime_byAlg];
                    
                    clear evnt_start LL_above_th evnt_indx evnt_lngth row1...
                        col1 events ev1 ev2 evnt_end ev_end n_to_combine  ...
                        SpMinMax SpAmp Abn0 Sp0 spikecount abnormalcount ...
                        seizurecount TotSPTime_byAlg TotABTime_byAlg ...
                        TotSZTime_byAlg sp short dse4_f dse5
                    
                    % events ev1 ev2
         
        else 
            seizures = [];
            abnormal  = [];
            spikes  = [];
            counts = [0 0 0];
            RealTotTime = [0 0 0];
            
    end
    
    total_seizures = [total_seizures; seizures];
    total_abnormal = [total_abnormal; abnormal];
    total_spikes = [total_spikes; spikes];

    
            % fill in structure with signal analysis
    
                % SigSumm(k,col).File = name;
%                 SigSumm(k,col).channel = col; 
%                 SigSumm(k,col).TotRecordLength = length(EEG)/f;
               % SigSumm(k,col).BulkSigScore = Thresholded(k);
%                 SigSumm(k,col).BaseStats = BaseStats(:,col);
%                 SigSumm(k,col).threshold = SF;
%                 SigSumm(k,col).EventCounts = counts;
%                 SigSumm(k,col).Spikes = spikes;
%                 SigSumm(k,col).Seizure = seizures;
%                 SigSumm(k,col).Abnormal = abnormal; 
%                 SigSumm(k,col).EventTimeTotals = RealTotTime;
                
                
                
        clear spHitssp spHitsabn spHitssz abnHitssp abnHitsabn ...
            abnHitssz szHitssp szHitsabn szHitssz spFNeg abnFNeg ...
            szFNeg SpikeFNeg AbnFNeg SeizureFNeg spFPos abnFPos ...
            szFPos SpikeFNeg AbnFNeg SeizureFNeg SpikeFPos ...
            AbnFPos SeizureFPos FNEGs FPOSs n_n spnc mdnc sznc spc mdc ...
            szc RealVid VidFNeg VidFPos_SP VidFPos_Abn VidFPos_SZ N_N ...
            ConfMat ConfMat2 SF2Cut SMinMax
        
    end

        clear SpEye AbnEye SzEye ByEye

end  

end

%% Detection pass 2

reset(ds);
ds.SelectedVariableNames = {'x8', 'x6', 'x31', 'x11', 'Time'}; % time and MSCEEG2 cage 8
ds.ReadSize = 100000;

EEG_waste = read(ds);

disp('seizure detection: pass 2')

while hasdata(ds)
    
    EEG_1 = table2array(read(ds));
    if hasdata(ds) EEG_2 = table2array(read(ds)); else EEG_2 = []; end;
    EEG = [EEG_1; EEG_2];
    progress = 100*(EEG(1,5)/time_end_8);
    elapsedTime = toc(timerVal);
    fprintf('Progress %0.1f %% ', progress)
    fprintf('in %0.1f min \r', elapsedTime/60)
    
for k = 1
%%
  
  
%     load (files(k).name);
%         input_name = files(k).name;            
%         [~, name, ~] = fileparts (input_name); 
        %EEG = clean;  
            fr=length(EEG)/f;
            t=EEG(:,5);             
            %t=t(:); 
    %%
    for col = 1:CH  % run all channels
                    
        [C,L]  = wavedec(EEG(:,col),DC,'db4');
          % perform wavelet decomposition to the indicated 
                    % decomposition level (DC) using the db4 wavelet
            L2 = [0; L];
             Decomp = C(1:L(1));
        %Decomp = C((sum(L2(1:ww-1))+1):L(ww));
            [rows, ~] = size (Decomp);
            % Set Decomp equal to the approximation at the DC level
       %%  
            % find distance between datapoints n and n-1
        Dist = zeros(rows,1);
            for r = 2:rows-w
                Dist(r-1) = abs(Decomp(r)-Decomp(r-1));
            end
            %%
            % sum distance between datapoints across window
            LLength = zeros(rows,1);
            LLStart = sum(Dist(1:w));
            LLength(1) = LLStart;
                for r = 2:rows-w
                    LL_a = LLength(r-1)-Dist(r-1);
                    LLength(r) = LL_a+Dist(r+w);
                end
            LLAve = LLength-BaseStats(1,col);
        
            clear rows LLAveSD
            [rows, ~] = size(Dist);
      
            LLAveSD = zeros(rows,1);
            LLAveSD = LLAve/(SF*BaseStats(2,col)); 
         
        % create binary vector of thresholded data
            [rows, columns] = size(LLAveSD);
                LLTh = zeros(rows,1);
                for r = 1:rows
                    if LLAveSD(r) >= 1
                        LLTh(r) = 1;
                    else, LLTh(r) = 0;
                    end
                end
                
        % erode (open) data to remove noise - half-window size = 1.
            LLTh2 = zeros(size(LLTh));
            nhood = 1;
            GG = length(LLTh);
                for r = 1
                    if sum(LLTh(1:(r+nhood))) == length(LLTh(1:(r+nhood)))
                        LLTh2(r) = 1;    
                    else, LLTh2(r) = 0;
                    end
                end
                for r = 2:GG-1
                    if sum(LLTh(r-nhood:r+nhood)) == nhood*2+1
                        LLTh2(r)=1;      
                    else, LLTh2(r) = 0;
                    end
                end
                for r = GG-1:GG
                    if sum(LLTh(r-nhood:GG)) == length(LLTh((r-nhood:GG)))
                        LLTh2(r) = 1;         
                    else, LLTh2(r) = 0;
                    end
                end
            
        % Dilate data to restore eroded points - half-window size = 1
            LLTh3 = zeros(size(LLTh));
                for r = 1
                    if sum(LLTh2(1:(r+nhood))) >= 1
                        LLTh3(r) = 1;      
                    else, LLTh3(r) = 0;
                    end
                end
                for r = 2:GG-1
                    if sum(LLTh2(r-nhood:r+nhood)) >=1
                        LLTh3(r)=1;          
                    else, LLTh3(r) = 0;
                    end
                end
                for r = GG-1:GG
                    if sum(LLTh2(r-nhood:GG)) >= 1
                        LLTh3(r) = 1; 
                    else, LLTh3(r) = 0;
                    end
                end
        
        % Dyadic upsampling at odd indices must be performed once for each
        % decomposition level.
            dse0 = LLTh3;
            dse1 = dyadup(dse0,1); %dyadic upsampling, odd indices (1)
            dse2 = dyadup(dse1,1); %dyadic upsampling, odd indices (2)
            dse3 = dyadup(dse2,1); %dyadic upsampling, odd indices (3)
            dse4 = dyadup(dse3,1); %dyadic upsampling, odd indices (4)
            dse5 = dyadup(dse4,1); 

           % if DC = 4,
           %     dse_U = dse4;
           % elseif DC = 5
           %     dse_U = dse5;
           % end
            
        % close events - half-window = half * decomp level squared to close
        % the effect of dyadic upsampling
            GG = length(dse4);
            nhood = 0.5*dc; 
            dse4_c = zeros(GG,2);
                for r = 1:nhood
                    if sum(dse4(1:(r+nhood))) ~= 0
                        dse4_c(r) = 1;
                        
                    else, dse4_c(r) = 0;
                    end
                end
                for r = nhood+1:GG-nhood
                    if sum(dse4(r-nhood:r+nhood)) ~=0
                        dse4_c(r)=1;
                    else, dse4_c(r) = 0;
                    end
                end
                for r = GG-nhood+1:GG
                    if sum(dse4(r-nhood:GG)) ~=0
                        dse4_c(r)=1;
                        
                    else, dse4_c(r) = 0;
                    end
                end

        % Erode the ends of events
                dse4_f = zeros(size(dse4_c));
                for r = 1:nhood
                    if sum(dse4_c(1:(r+nhood))) == ...
                            length(dse4_c(1:(r+nhood)))
                        dse4_f(r) = 1;
                    else ,dse4_f(r) = 0;
                    end
                end
                for r = nhood+1:GG-nhood
                    if sum(dse4_c(r-nhood:r+nhood)) == nhood*2+1
                        dse4_f(r)=1;
                    else ,dse4_f(r) = 0;
                    end
                end
                for r = GG-nhood+1:GG
                    if sum(dse4_c(r-nhood:GG)) == ...
                            length(dse4_c((r-nhood):GG))
                        dse4_f(r) = 1;
                    else ,dse4_f(r) = 0;
                    end
                end
        
            LL_above_th_row = find(dse4_f==1); % edits to convert row to time 2/2 datastore
            LL_above_th = EEG(LL_above_th_row,5);    
                
                % pick out line lengths greater than 
                % BaseMedian + SF*StDev. return as 
                % row index (time) of event start
                       
            clear LLAve LLAveSD LLTh LLTh2 LLTh3 nhood GG r dse0 dse1 ...
                dse2 dsde3 dse4 dse4_c 
 %%               
    if isempty(LL_above_th)~=1
        
        %% Event classification    
        
        % RJK event times
        
        % organize all threshold crossings with start, end, length
        all_above_th = zeros((length(LL_above_th)-1), 3);           
        all_above_th(:,1) = LL_above_th(1:(length(LL_above_th)-1),1);
        all_above_th(:,2) = LL_above_th(2:end,1);
        all_above_th(:,3) = all_above_th(:,2) - all_above_th(:,1);
        
        % identify and organize events as those with > 750ms intervals
        % between
        if size(all_above_th,1) == 0, events = [0 0 0];
        elseif size(all_above_th,1) == 1, events = all_above_th;
        else events = zeros(length(find(all_above_th(:,3) > .75)),3); ...   
            events(:,2) = all_above_th(find(all_above_th(:,3) > 0.75),1); ...
            events(:,1) = all_above_th((find(all_above_th(:,3) > 0.75)),2); ...
            events(:,1) = circshift(events(:,1),1); ...
            events(1,1) = all_above_th(1,1);
        end
        
        % if size(events,1) >1 events(2:end,1) = events(1:(size(events,1)-1),1); ...
        %    elseif size(events,1) ~= 0 events(1,1) = all_above_th(1,1);, else events = [0 0 0];,end;
        %if length(events) ~= 0 events(:,3) = events(:,2) - events (:,1);, else events = [];,end;
        
            % identify hits that are greater than 1 from neighboring event
%                evnt_start = zeros(length(LL_above_th),1);  
%                     for r = 2:length(LL_above_th)
%                         evnt_start(r,1) = LL_above_th(r) - LL_above_th(r-1);
%                     end
%                 evnt_indx = find(evnt_start > 0.001); % edited for actual time vs rows
%        
%         % Find event lengths
%                 evnt_lngth = zeros(length(evnt_indx),1);
%                     for r = 1:length(evnt_indx)-1
%                         evnt_lngth(r,1) = LL_above_th(evnt_indx(r+1)-1) ...
%                             - LL_above_th(evnt_indx(r));
%                     end
%                 evnt_lngth(length(evnt_indx),1) = LL_above_th(end,1)...
%                     - LL_above_th(evnt_indx(end,1));
% 
%                 [row1, col1] = size(evnt_lngth);
% 
%         % Find event ends
%                 evnt_end = zeros(size(evnt_indx));
%                     for r = 1:length(evnt_indx-1)
%                         evnt_end(r,1) = LL_above_th(evnt_indx(r))+ evnt_lngth(r);
%                     end
%                 ev_end = [LL_above_th(evnt_indx) evnt_end];
%               %  ev_end(1,3) = 400;
%                 ev_end(1,3) = f;
%                     for r = 2:length(evnt_indx)
%                             ev_end(r,3) = ev_end((r),1)-ev_end((r-1),2);
%                     end
%     
%         % Identify events less than 750 ms from neighbor.
%         % Events closer than 750 ms to neighbor will be combined with neighbor.
%                 ev_end(:,4) = 0;
%                     for r=1:length(evnt_indx)
%                         if ev_end(r,3) <= .75*f
%                            ev_end(r,4) = 0;
%                         else ,ev_end(r,4) = 1;
%                         end
%                     end
%                 n_to_combine = find(ev_end(:,4) == 1); 
%                     % index of events closer than 250 ms
%     
%         % recreate event list after combining neighbors 
%                 events = zeros(length(n_to_combine),2);
%                     for r = 1:length(n_to_combine)-1
%                         events(r,1) = ev_end(n_to_combine(r),1);
%                         events(r,2) = ev_end((n_to_combine(r+1)-1),2);
%                     end
%                 events(length(n_to_combine),1) = ...
%                     ev_end(n_to_combine(end),1);
%                 events(length(n_to_combine),2) = ev_end(end,2);
    
        % Sort events based on event lengths
            % Short = events < 3 sec, includes spikes and abnormal events; 
            % Seizures = events > 3 sec in length 
                events(:,3) = events(:,2)-events(:,1);
                [rows, columns] = size(events);
                    for r = 1:rows
                        %if events(r,3) < 5*f
                         if events(r,3) < 3 % adjusted from 1 (sec?) to 3
                            events(r,4) = 1;
                        else ,events(r,4) = 2;
                        end
                    end
                ev1 = find(events(:,4)==1);
                ev2 = find(events(:,4)==2);

                    
                    if isempty(ev2) ~= 1
                        seizures = events(ev2,1:2); % removed division by f
                    else ,seizures = [];
                    end

        % Sort Short events based on event amplitude
            % Spikes = events > 250 uV in amplitude
            % Abnormal = events < 250 uV in amplitude
                if isempty(ev1) ~= 1 
                    short = events(ev1,1:2);
                    [rows, ~] = size(short);
                    SpMin = zeros(rows,2);
                    
                   EEGev = zeros(length(dse4_f),2);
                   EEGev(1:length(EEG),2) = EEG(:,5); % added time column to EEGev
                   EEGev(1:length(EEG),1) = EEG(:,col);
                    for r = 1:rows
                        sp = EEGev(find(short(r,1)==EEGev(:,2)):find(short(r,2)==EEGev(:,2)),1);
                        % sp = sp';
                        if sp ~=0 SpMinMax(r,:) = [min(sp) max(sp)];, else SpMinMax = zeros(r,2); end;
                    end
        
                SpAmp = SpMinMax(:,2)-SpMinMax(:,1);
                    Abn0 = find(SpAmp <.250);
                        abnormal = short(Abn0,1:2);
                    Sp0 = find(SpAmp>=.250);
                        spikes = short(Sp0,1:2);    
                else
                        spikes = [];
                        abnormal = [];
                end

                spikecount = length(spikes);
                abnormalcount = length(abnormal);
                seizurecount = length(seizures);
                    counts = [spikecount seizurecount abnormalcount ];
          
     % find the actual abnormal time totals
    
                if isempty(spikes) ~=1
                    TotSPTime_byAlg = sum(spikes(:,2)-spikes(:,1));
                else ,TotSPTime_byAlg = 0;
                end
                if isempty(abnormal)~=1
                    TotABTime_byAlg = sum(abnormal(:,2)-abnormal(:,1));
                else ,TotABTime_byAlg = 0;
                end
                if isempty(seizures) ~= 1
                    TotSZTime_byAlg = sum(seizures(:,2)-seizures(:,1));
                else ,TotSZTime_byAlg = 0;
                end
                    RealTotTime = [TotSPTime_byAlg  ...
                        TotSZTime_byAlg TotABTime_byAlg];
                    
                    clear evnt_start LL_above_th evnt_indx evnt_lngth row1...
                        col1 events ev1 ev2 evnt_end ev_end n_to_combine  ...
                        SpMinMax SpAmp Abn0 Sp0 spikecount abnormalcount ...
                        seizurecount TotSPTime_byAlg TotABTime_byAlg ...
                        TotSZTime_byAlg sp short dse4_f dse5
                    
                    % events ev1 ev2
         
        else 
            seizures = [];
            abnormal  = [];
            spikes  = [];
            counts = [0 0 0];
            RealTotTime = [0 0 0];
            
    end
    
    total_seizures = [total_seizures; seizures];
    total_abnormal = [total_abnormal; abnormal];
    total_spikes = [total_spikes; spikes];

    
            % fill in structure with signal analysis
    
                % SigSumm(k,col).File = name;
%                 SigSumm(k,col).channel = col; 
%                 SigSumm(k,col).TotRecordLength = length(EEG)/f;
               % SigSumm(k,col).BulkSigScore = Thresholded(k);
%                 SigSumm(k,col).BaseStats = BaseStats(:,col);
%                 SigSumm(k,col).threshold = SF;
%                 SigSumm(k,col).EventCounts = counts;
%                 SigSumm(k,col).Spikes = spikes;
%                 SigSumm(k,col).Seizure = seizures;
%                 SigSumm(k,col).Abnormal = abnormal; 
%                 SigSumm(k,col).EventTimeTotals = RealTotTime;
                
                
                
        clear spHitssp spHitsabn spHitssz abnHitssp abnHitsabn ...
            abnHitssz szHitssp szHitsabn szHitssz spFNeg abnFNeg ...
            szFNeg SpikeFNeg AbnFNeg SeizureFNeg spFPos abnFPos ...
            szFPos SpikeFNeg AbnFNeg SeizureFNeg SpikeFPos ...
            AbnFPos SeizureFPos FNEGs FPOSs n_n spnc mdnc sznc spc mdc ...
            szc RealVid VidFNeg VidFPos_SP VidFPos_Abn VidFPos_SZ N_N ...
            ConfMat ConfMat2 SF2Cut SMinMax
        
    end

        clear SpEye AbnEye SzEye ByEye

end  

end

elapsedTime = toc;
fprintf('Total time %0.1f min \r', elapsedTime/60)
  
%  save ('A1_SigSumm_Ch1-8.mat', 'SigSumm')
