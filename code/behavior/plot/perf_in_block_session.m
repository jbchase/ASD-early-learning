function perf = perf_in_block_session(results, protocol, blockLength)

% calculate performance in trial blocks but on a session-basis

% AB: calculate 3 sessions
% CD: calculate 3 sessions
% DC: calculate 6 sessions

if strcmp(protocol, 'AB') || strcmp(protocol,'CD') || strcmp(protocol, 'AB-CD')
    nSessions = 3;
elseif strcmp(protocol, 'DC') || strcmp(protocol, 'AB-DC')
    nSessions = 6;
end

perf = nan(1,20, nSessions);

for ss = 1:nSessions
    nTrials = size(results{ss},1);
    for qq = 1:20
        startTrial = blockLength*(qq-1) + 1;
        endTrial = blockLength*qq;
        if endTrial < nTrials && startTrial < nTrials
            perf(1,qq,ss) = sum(~isnan(results{ss}.reward(startTrial:endTrial)))/blockLength;
        elseif endTrial > nTrials && startTrial < nTrials
            perf(1,qq,ss) = sum(~isnan(results{ss}.reward(startTrial:end)))/(nTrials-startTrial+1);
        end
    end
end
