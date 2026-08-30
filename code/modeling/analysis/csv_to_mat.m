function X = csv_to_mat(csvFile, matFile, scheduleType)
%CSV_TO_MAT  Convert a raw behavior CSV export to a processed .mat table.
%
%   X = CSV_TO_MAT(csvFile, matFile, scheduleType) reads the raw behavior
%   CSV exported from the operant-box system, recodes it into the analysis
%   table format used by the modeling pipeline, and saves it to matFile.
%
%   scheduleType is 'deterministic' (every correct choice rewarded) or
%   'probabilistic' (10-20% of correct choices unrewarded).
%
%   The processed table X contains the columns used by the pipeline:
%     Animal_ID - numeric subject ID
%     genotype  - 0 = WT, 1 = Het, 2 = KO
%     schedule  - odor cue: 1 = A, 2 = B (recoded; values >= 3 are removed)
%     action    - chosen port: 1 = left, 2 = right
%     reward1   - 1 = rewarded, 0 = unrewarded
%
%   The processed .mat files in data/processed are the canonical inputs of
%   the modeling pipeline; this function documents how they were derived
%   from the raw CSVs in data/raw.
%
%   See also FIT_MODELS.

%   Author: Jing-Jing Li (jl3676@berkeley.edu)

X = readtable(csvFile);

% Female exports use "results"/"schedules" column names
if ismember('results', X.Properties.VariableNames)
    X = renamevars(X, ["results", "schedules"], ["result", "schedule"]);
end

switch scheduleType
    case 'deterministic'
        % Keep completed odor-learning trials; recode result to reward1
        X(X.result > 2 | X.result <= 1, :) = [];
        X.result(X.result == 2) = 0;
        X.result(X.result == 1.2) = 2;
        X.result(X.result == 1.3) = 3;
        X = renamevars(X, ["result", "geno"], ["reward1", "genotype"]);

        % Numeric subject IDs (from the 4th-6th characters of animal_ids)
        X.Animal_ID = cellfun(@(s) str2double(s(4:6)), X.animal_ids);
        X.Animal_ID(cellfun(@(s) contains(s, 'TSC2'), X.animal_ids)) = 6;

        % Recode the correct port into the chosen action
        X.action(X.reward1 > 0) = mod(X.schedule(X.reward1 > 0) - 1, 2) + 1;
        X.action(X.reward1 == 0) = mod(3 - X.schedule(X.reward1 == 0) - 1, 2) + 1;

    case 'probabilistic'
        X = renamevars(X, ["geno"], ["genotype"]);
        X.schedule = mod(X.schedule - 1, 2) + 1;
        X.action(X.result == 1.1 | X.result == 1.3 | X.result == 1.2) = X.schedule(X.result == 1.1 | X.result == 1.3 | X.result == 1.2);
        X.action(X.result == 2) = 3 - X.schedule(X.result == 2);
        X.reward1(X.portsides > 0 & (X.result == 1.3 | X.result == 1.2)) = 1;
        X.reward1(X.portsides == 0 | X.result == 2) = 0;
        X.reward1 = X.reward1 > 0;

    otherwise
        error('csv_to_mat:unknownSchedule', 'Unknown schedule type: %s', scheduleType);
end

save(matFile, 'X');
end
