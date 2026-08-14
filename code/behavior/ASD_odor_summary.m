function ASD_odor_summary(dataIndex, strain, savefigpath, savedatapath)

%% go over every session to plot performance in blocks for
% 1. first two sessions of AB
% 2. first two sessions of CD
% 3. DC reverse

%% number of AB/CD/DC trials experienced based on training session
    % training session: 3 AB + 3 AB-CD + 6 AB-(CD) - DC sessions
    % ignore the retrain sessions for now

nFiles = size(dataIndex,1);
Subjects = unique(dataIndex.Animal);
nSubjects = length(unique(dataIndex.Animal));

genotype = cell(nSubjects,1);

% number of trials
nAB_trials = nan(nSubjects, 12);
nCD_trials = nan(nSubjects, 4);
nDC_trials = nan(nSubjects, 6);

% trial timing
timing = 1:8*60; % in minutes (8 hour maximum)
perf_timing = 1:30:8*60;
Trials_time = nan(nFiles, length(timing));
performance_time = nan(nFiles, length(perf_timing));

% performance
blockLength = 100; % check performance in blocks of 100 trials
perf_AB_block = nan(nSubjects, 50);
perf_CD_block = nan(nSubjects, 50);
perf_DC_block = nan(nSubjects, 50);
perf_AB_CD_block = nan(nSubjects, 4);  % AB readiness performance in AB-CD session
perf_AB_DC_block = nan(nSubjects, 4); % AB readiness performance in AB-DC session

perf_AB_quantile = nan(nSubjects,4, 3); % dim3: first session, second session
perf_CD_quantile = nan(nSubjects,4, 3); % dim3: session
perf_DC_quantile = nan(nSubjects,4, 6); % dim4: session

perf_AB_3h = nan(nSubjects, 20, 3); % get trials in the first 3h

perf_AB_block_session = nan(nSubjects,20, 3); % dim3: first session, second session
perf_CD_block_session = nan(nSubjects,20, 3); % dim3: session
perf_DC_block_session = nan(nSubjects,20, 6); % dim4: session
perf_AB_CD_block_session = nan(nSubjects, 20, 3);  % AB readiness performance in AB-CD session
perf_AB_DC_block_session = nan(nSubjects, 20, 6); % AB readiness performance in AB-DC session
perf_AB_block_learning= nan(nSubjects, 20, 3);
perf_AB_running = nan(nSubjects,2000,3);
perf_AB_running_smoothed = nan(nSubjects,2000,3);
perf_DC_session = nan(nSubjects,6);

%% trial number and performance
for ii = 1:nSubjects

    genotype(ii) = unique(dataIndex.Genotype(strcmp(dataIndex.Animal, Subjects{ii})));
    %% look for first 3 AB sessions
    subDataIndex = dataIndex(strcmp(dataIndex.Animal, Subjects{ii}),:);
    % load the first three sessions
    protocol = 'AB';

    nSessions = 3;
    results = cell(nSessions,1);
    for ss=1:nSessions
        if ss<= size(subDataIndex,1)
            csvFile = fullfile(subDataIndex.BehPath{ss}, ...
                subDataIndex.BehCSV{ss});
            if exist(csvFile)
                results{ss} = readtable(csvFile);
                nAB_trials(ii, ss) = size(results{ss},1);
            else
                results{ss} = NaN;
                nAB_trials(ii,ss) = NaN;
            end
        else
            results{ss} = NaN;
            nAB_trials(ii,ss) = NaN;
        end
    end
    % calculate performance in quantile
    %perf_AB_quantile(ii,:,:) = perf_in_quantile(results, protocol);
    
    % calculate performance in 100-trial block on a session basis
    perf_AB_block_session(ii,:,:) = perf_in_block_session(results, protocol, blockLength);
    
    % calculate performance in 100-trial block within the first 3 hours
    results_3h = {};
    for rii = 1:length(results)
        tempResults = results{rii};
        if size(tempResults,1)>1
            results_3h{rii} = tempResults((results{rii}.center_in-results{rii}.center_in(1))<3*60*60,:);
        else
            results_3h{rii} = [NaN];
        end
    end
    perf_AB_3h(ii,:,:) = perf_in_block_session(results_3h, protocol, blockLength);
    % calculate performance in blocks
    %perf_AB_block(ii,:) = perf_in_block(results, protocol, blockLength);
    % examine the learning rate of AB sessions with 100 block length
    perf_AB_running(ii,:,:) = perf_in_running_session(results,protocol,blockLength);
    
%     
%     order = 2;       % polynomial order
%     framelen = 21;   % must be odd
%     for ss = 1:3
%         tempSig = perf_AB_running(ii,:,ss);
%         tempSig = tempSig(~isnan(tempSig));
%         perf_AB_running_smoothed(ii,1:length(tempSig),ss) = sgolayfilt(tempSig, order, framelen);
%     end
    
    %% for CD odors
    animalMask = strcmp(dataIndex.Animal, Subjects{ii});
    stageMask = cellfun(@(x) isequal(x, [1;2;3;4]), dataIndex.OdorPresented);
    subDataIndex = dataIndex(animalMask & stageMask,:);

    protocol = 'CD';
    nSessions = 3;
    results = cell(nSessions,1);
    results_AB = cell(nSessions,1);
    for ss=1:nSessions
        if ss<= size(subDataIndex,1)
            csvFile = fullfile(subDataIndex.BehPath{ss}, ...
                subDataIndex.BehCSV{ss});

            results{ss} = readtable(csvFile);       

            % find the switch point
            switchTrial = find(results{ss}.schedule == 3 | results{ss}.schedule == 4, 1);
            results_AB{ss} = results{ss}(1:switchTrial-1,:);
            results{ss} = results{ss}(switchTrial:end,:);
            
            nAB_trials(ii, ss+3) = switchTrial-1;
            nCD_trials(ii, ss) = size(results{ss},1);
        else
            results{ss} = NaN;
            results_AB{ss} = NaN;
        end
    end
    % calculate performance in quantile
    perf_CD_quantile(ii,:,:) = perf_in_quantile(results, protocol);
    
     % calculate performance in quantile on a session basis
    perf_CD_block_session(ii,:,:) = perf_in_block_session(results, protocol, blockLength);

    perf_AB_CD_block_session(ii,:,:) = perf_in_block_session(results_AB, 'AB-CD', blockLength);
    % calculate performance in blocks
    %perf_CD_block(ii,:) = perf_in_block(results, protocol,blockLength);

    %% for DC performance

    animalMask = strcmp(dataIndex.Animal, Subjects{ii});
    stageMask = cellfun(@(x) isequal(x, [1;2;5;6]), dataIndex.OdorPresented) | cellfun(@(x) isequal(x, [1;2;3;4;5;6]), dataIndex.OdorPresented);
    subDataIndex = dataIndex(animalMask & stageMask,:);

    if ~isempty(subDataIndex)
        protocol = 'DC';
        nSessions = 6;
        results = cell(nSessions,1);
        results_AB = cell(nSessions, 1);
        for ss =1 :nSessions
            if ss<= size(subDataIndex,1)
                csvFile = fullfile(subDataIndex.BehPath{ss}, ...
                    subDataIndex.BehCSV{ss});

                results{ss} = readtable(csvFile);

                % find the switch point
                switchCD =  find(results{ss}.schedule == 3 | results{ss}.schedule == 4, 1);
                switchTrial = find(results{ss}.schedule == 5 | results{ss}.schedule == 6, 1);
                
                if switchCD > 0
                    nAB_trials(ii, ss+6) = switchCD-1;
                    nCD_trials(ii,ss+3) = switchTrial-switchCD;
                    results_AB{ss} = results{ss}(1:switchCD-1,:);
                else
                    nAB_trials(ii,ss+6) = switchTrial-1;
                    results_AB{ss} = results{ss}(1:switchTrial-1,:);
                end
                results{ss} = results{ss}(switchTrial:end,:);

                nDC_trials(ii, ss) = size(results{ss},1);
            else
                results{ss} = NaN;
                results{ss} = NaN;
            end
        end
        % calculate performance in quantile
        perf_DC_quantile(ii,:,:) = perf_in_quantile(results, protocol);
        
         % calculate performance in quantile on a session basis
        perf_DC_block_session(ii,:,:) = perf_in_block_session(results, protocol, blockLength);
        perf_AB_DC_block_session(ii,:,:) = perf_in_block_session(results_AB, 'AB-DC', blockLength);

        % calculate performance in blocks
        %perf_DC_block(ii,:) = perf_in_block(results, protocol,blockLength);
    end
end

%% calculate slope of smoothed performance curve
windowSize = 401;   % must be odd
dx = nan(size(perf_AB_running_smoothed));
halfWin = floor(windowSize/2);

for ss = 1:size(perf_AB_running_smoothed,1)
    for pp = 1:size(perf_AB_running_smoothed,3)
for i = 1+halfWin : size(perf_AB_running_smoothed,2)-halfWin
    y = perf_AB_running_smoothed(ss,i-halfWin : i+halfWin,pp);
    t = (1:windowSize)';         % local time index
    p = polyfit(t, y, 1);        % linear fit
    dx(ss,i,pp) = p(1);                % slope
end
    end
end

%% check DC performance based on how many AB/CD trials they have experienced
% check every session before starting DC reversal, including retrain
% session

nAB_trials_total = {};
nCD_trials_total = {};
nABinDC_trials_total = {};
for ii = 1:nSubjects
    %% look for first 3 AB sessions
    subDataIndex = dataIndex(strcmp(dataIndex.Animal, Subjects{ii}),:);
    tempTrial_AB = [];
    tempTrial_CD = [];
    tempTrial_ABinDC = [];
    reversed = 0; % check if the animal had experienced reversal already
    % used to identify AB retrain sessions in during reversal
    for ss = 1:length(subDataIndex.Protocol)
        csvFile = fullfile(subDataIndex.BehPath{ss}, ...
            subDataIndex.BehCSV{ss});
        if exist(csvFile)
            result = readtable(csvFile);
        if ~ strcmp('AB-CD-DC', subDataIndex.Protocol{ss})
                if reversed==0
                    Index_AB = find(result.schedule==3 | result.schedule==4,1);
                    if isempty(Index_AB)
                        tempTrial_AB = [tempTrial_AB, size(result,1)];
                    else
                        tempTrial_AB = [tempTrial_AB, Index_AB-1];
                        Index_CD = find(result.schedule==3|result.schedule==4, 1, 'last');
                        tempTrial_CD = [tempTrial_CD, Index_CD-Index_AB+1];
                    end
                else
                    Index_AB = find(result.schedule==5 | result.schedule==6,1);
                    if isempty(Index_AB)
                        tempTrial_ABinDC = [tempTrial_ABinDC, size(result,1)];
                    else
                        tempTrial_ABinDC = [tempTrial_ABinDC, Index_AB-1];
                    end
                end
                elseif strcmp('AB-CD-DC', subDataIndex.Protocol{ss})
                    reversed =1;
                    Index_AB = find(result.schedule==3 | result.schedule==4,1);
                    if isempty(Index_AB)
                    tempTrial_AB = [tempTrial_AB, size(result,1)];
                else
                    tempTrial_AB = [tempTrial_AB, Index_AB-1];
                    Index_CD = find(result.schedule==3|result.schedule==4, 1, 'last');
                    tempTrial_CD = [tempTrial_CD, Index_CD-Index_AB+1];
                end
            end
        else
            tempTrial_AB=NaN;
            tempTrial_CD = NaN;
            tempTrial_DC = NaN;
        end
    end
    nAB_trials_total{ii}=tempTrial_AB;
    nCD_trials_total{ii} = tempTrial_CD;
    nABinDC_trials_total{ii} = tempTrial_ABinDC;
end

%% trial timing and performance
for ss=1:nFiles
   csvFile = fullfile(dataIndex.BehPath{ss}, ...
                dataIndex.BehCSV{ss});
   if exist(csvFile)
   results = readtable(csvFile);
   for tt = timing
       startSeconds = (tt-1)*60;
       endSeconds = tt*60-1;
       sessionTiming = results.center_in-results.center_in(1);
       %if endSeconds <= sessionTiming(end)
       trialMask=(sessionTiming>=startSeconds & sessionTiming<=endSeconds);
       %elseif endSeconds > sessionTiming(end) && startSeconds <=sessionTiming(end)
       %     trialMask=(sessionTiming>=startSeconds);
       %else

       Trials_time(ss,tt) = sum(trialMask);
   end

   for pp = 1:length(perf_timing)
       startSeconds = (pp-1)*30*60;
       endSeconds = pp*30*60-1;
       sessionTiming = results.center_in-results.center_in(1);
       trialMask=(sessionTiming>=startSeconds & sessionTiming<=endSeconds);
       performance_time(ss,pp) = sum(~isnan(results.reward(trialMask)))/sum(trialMask);
   end
   else
       Trials_time(ss,:) = NaN;
       performance_time(ss,:) = NaN;
   end

end

%% make the plot

setup_figprop;

%% percentage of total trials performed by time
Trials_time_percentage = cumsum(Trials_time,2)./sum(Trials_time,2);
%plot_Trial_distrubution(Trials_time, performance_time, savefigpath);
% cap at 60 minues
%% number of trials performed

plot_nTrials(nAB_trials, nCD_trials, nDC_trials,genotype, savefigpath);

%% performance in block within the first 3 hours
plot_performance_block_session(perf_AB_3h, strain, genotype, savefigpath, 'AB-3h');

%% performance in block on session-basis
plot_performance_block_session(perf_AB_block_session, strain, genotype, savefigpath, 'AB');

%% save the data in csv format
data_s1 = squeeze(perf_AB_block_session(:,:,1));   % 20 x 20

T1 = array2table(data_s1);
T1.Genotype = genotype(:);                         % align as column
T1.Subjects
T1 = movevars(T1,'Genotype','Before',1);           % put genotype first
savecsvpath = fullfile(savedatapath, [strain, 'session1.csv']);
writetable(T1,savecsvpath);


% Session 2
data_s2 = squeeze(perf_AB_block_session(:,:,2));   % 20 x 20

T2 = array2table(data_s2);
T2.Genotype = genotype(:);
T2 = movevars(T2,'Genotype','Before',1);
savecsvpath = fullfile(savedatapath, [strain, 'session2.csv']);
writetable(T2,savecsvpath);

data_s1 = squeeze(perf_AB_3h(:,:,1));   % 20 x 20

T1 = array2table(data_s1);
T1.Genotype = genotype(:);                         % align as column
T1 = movevars(T1,'Genotype','Before',1);           % put genotype first
savecsvpath = fullfile(savedatapath, [strain, 'session1-3h.csv']);
writetable(T1,savecsvpath);


% Session 2
data_s2 = squeeze(perf_AB_3h(:,:,2));   % 20 x 20

T2 = array2table(data_s2);
T2.Genotype = genotype(:);
T2 = movevars(T2,'Genotype','Before',1);
savecsvpath = fullfile(savedatapath, [strain, 'session2-3h.csv']);
writetable(T2,savecsvpath);