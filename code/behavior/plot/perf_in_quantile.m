function perf = perf_in_quantile(results, protocol)

% AB: calculate 3 sessions
% CD: calculate 3 sessions
% DC: calculate 6 sessions

if strcmp(protocol, 'AB') || strcmp(protocol,'CD')
    nSessions = 3;
elseif strcmp(protocol, 'DC')
    nSessions = 6;
end

perf = nan(1,4, nSessions);
nTrials = cell(nSessions,1);
for ss = 1:nSessions
    nTrials{ss} = size(results{ss},1);

    for qq = 1:4
        startTrial = floor(nTrials{ss}/4)*(qq-1) + 1;
        endTrial = floor(nTrials{ss}/4)*qq;
        if size(results{ss},1)>1
            perf(1,qq, ss) = sum(~isnan(results{ss}.reward(startTrial:endTrial)))/(endTrial-startTrial+1);
        end
    end
end
