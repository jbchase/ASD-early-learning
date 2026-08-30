function R = validate_models(id, opts)
%VALIDATE_MODELS  Simulate the winning model and compare against behavior.
%
%   R = VALIDATE_MODELS(id) simulates choice behavior for each animal of
%   cohort ID 1000 times using the MAP parameters of the winning model
%   (a0b1s_hybrid), using the latent state occupancy trajectories inferred
%   from the animal's data, and compares simulated and observed learning
%   curves at the individual and group level (manuscript Fig 2C-E,
%   Supplementary Figs 5B,C,E,F and 7B,C,E,F). Per-trial p(engaged)
%   trajectories are saved to the cohort's p(engaged) file for cohorts
%   where the manuscript reports engagement analyses.
%
%   VALIDATE_MODELS(id, opts) allows optional fields:
%     nIters         - simulations per animal (default: 1000)
%     subjectIndices - validate only these subject indices (default: all)
%     fitFile        - load fits from this file instead of the cohort's
%                      archived results file
%     seed           - if nonempty, call rng(seed) before simulating and run
%                      the simulation loop serially, making output
%                      reproducible
%     makePlots      - make validation figures (default: true)
%     saveResults    - save the p(engaged) file (default: true)
%
%   The returned struct R contains:
%     allLatent          - 1 x nSubjects cell array of p(engaged) trajectories
%     allBehavQuartiles  - 4 x nSubjects observed quartile accuracies
%     allSimQuartiles    - 4 x nSubjects simulated quartile accuracies
%     groups             - group index per subject (1 = Het, 2 = WT, 3 = KO)
%
%   See also COHORTS, FIT_MODELS, COMPARE_PARAMS.

%   Author: Jing-Jing Li (jl3676@berkeley.edu)

arguments
    id (1, :) char
    opts.nIters (1, 1) double = 1000
    opts.subjectIndices = []
    opts.fitFile = ''
    opts.seed = []
    opts.makePlots (1, 1) logical = true
    opts.saveResults (1, 1) logical = true
end

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'models', 'likelihood'));
addpath(fullfile(repoRoot, 'models', 'generative'));

C = cohorts(id);
model_name = 'a0b1s_hybrid';
windowSize = 20;

if isempty(opts.fitFile)
    opts.fitFile = fullfile(repoRoot, 'results', C.fitFile);
end
m = load(opts.fitFile);
X = m.X; subjects = m.subjects; Ms = m.Ms; All_Params = m.All_Params; All_fits = m.All_fits;
model_ind = strcmp(cellfun(@(x) x.name, Ms(:), 'UniformOutput', false), model_name);
AICs = squeeze(All_fits(:, 3, :));
if ~isempty(opts.subjectIndices)
    subjects = subjects(opts.subjectIndices);
    AICs = AICs(opts.subjectIndices, :);
    All_Params{model_ind} = All_Params{model_ind}(opts.subjectIndices, :);
end

%% Assign animal groups (1 = Het, 2 = WT, 3 = KO)
switch C.groupScheme
    case 'two'
        groups = zeros(length(subjects), 1);
        for k = 1:length(subjects)
            T = find((X.Animal_ID == subjects(k)) & X.schedule < 3);
            if X.genotype(T(1)) == 1
                groups(k) = 1;
            else
                groups(k) = 2;
            end
        end
    case 'three'
        groups = zeros(length(subjects), 1);
        for k = 1:length(subjects)
            T = find((X.Animal_ID == subjects(k)) & X.schedule < 3);
            if X.genotype(T(1)) == 2
                groups(k) = 3;      % KO
            elseif X.genotype(T(1)) == 1
                groups(k) = 1;      % Het
            else
                groups(k) = 2;      % WT
            end
        end
end
groupIds = unique(groups)';

%% Per-subject behavior, latent states, and simulations
if ~isempty(opts.seed)
    rng(opts.seed);
end
allLatent = cell(1, length(subjects));
allBehavQuartiles = zeros(4, length(subjects));
allSimQuartiles = zeros(4, length(subjects));

if opts.makePlots
    fig = figure('Position', [300, 300, 1500, 800], 'Visible', 'off');
end

for k = 1:length(subjects)
    s = subjects(k);
    this_data = get_subject_data(X, s);

    quartiles = round(quantile(1:size(this_data, 1), 3));
    allBehavQuartiles(:, k) = [nanmean(this_data(1:quartiles(1), 3)), ...
        nanmean(this_data(quartiles(1):quartiles(2), 3)), ...
        nanmean(this_data(quartiles(2):quartiles(3), 3)), ...
        nanmean(this_data(quartiles(3):end, 3))];

    % Simulate the model using the inferred latent state trajectory
    theta = All_Params{model_ind}(k, :)';
    latent = feval([model_name, '_latent'], theta, this_data);
    allLatent{k} = latent;
    sim_rewards = zeros(size(this_data, 1), 1);
    if isempty(opts.seed)
        % Production setting: parallel simulation loop
        parfor it = 1:opts.nIters
            sim_rewards = sim_rewards + simulate_once(model_name, theta, this_data, latent, C.pReward);
        end
    else
        % Test setting: serial, seeded, reproducible
        for it = 1:opts.nIters
            sim_rewards = sim_rewards + simulate_once(model_name, theta, this_data, latent, C.pReward);
        end
    end
    sim_rewards = sim_rewards / opts.nIters;
    allSimQuartiles(:, k) = [mean(sim_rewards(1:quartiles(1))), ...
        mean(sim_rewards(quartiles(1):quartiles(2))), ...
        mean(sim_rewards(quartiles(2):quartiles(3))), ...
        mean(sim_rewards(quartiles(3):end))];

    if opts.makePlots
        figure(fig);
        subplot(C.indivSubplots(1), C.indivSubplots(2), k)
        plot(movmean(this_data(:, 3), windowSize, "omitnan"), 'Color', 'blue', 'LineWidth', 2)
        ylabel('Accuracy')
        ylim([0 1])
        if groups(k) == 3, genotype = 'KO';
        elseif groups(k) == 1, genotype = 'HT';
        else, genotype = 'WT'; end
        title([genotype num2str(s) ' (' num2str(round(exp(-(AICs(k, model_ind) - 2*length(theta)) / (2*size(this_data, 1))), 2)) ')'])
        hold on
        plot(movmean(latent, windowSize/2), 'Color', 'black')
        hold on
        plot(movmean(sim_rewards, windowSize/2, 'omitnan'), 'Color', 'red', 'LineWidth', 2)
        xlim([0 size(sim_rewards, 1)])
    end
end

%% Save p(engaged) trajectories
if opts.saveResults && C.pEngaged
    save(fullfile(repoRoot, 'results', C.pEngagedFile), 'allLatent');
end

R = struct('cohort', id, 'allLatent', {allLatent}, 'groups', groups, ...
    'allBehavQuartiles', allBehavQuartiles, 'allSimQuartiles', allSimQuartiles);

%% Figures
if opts.makePlots
    figure(fig);
    sgtitle([id ' Sess1 ' model_name], 'Interpreter', 'none')
    plotDir = fullfile(repoRoot, 'plots');
    if ~exist(plotDir, 'dir'), mkdir(plotDir); end
    saveas(fig, fullfile(plotDir, [id '_validate_individual.svg']))
    close(fig)

    % Group validation: observed vs simulated quartile accuracies
    groupColors = {'#D95319', 'blue', 'black'};
    fig = figure('Position', [300, 300, 600, 400], 'Visible', 'off');
    plot([1, 4], [.5 .5], '--', 'Color', 'black')
    hold on
    plot([1, 4], [.7 .7], '--', 'Color', 'black')
    hold on
    legEntries = {};
    for g = groupIds
        gd = allBehavQuartiles(:, groups == g);
        gs = allSimQuartiles(:, groups == g);
        meanData = mean(gd, 2);
        semData = std(gd, 0, 2) ./ sqrt(size(gd, 2));
        meanSimData = mean(gs, 2);
        errorbar(1:4, meanData, semData, '-o', 'Color', groupColors{g}, 'LineWidth', 2, 'MarkerFaceColor', groupColors{g})
        hold on
        plot(1:4, meanSimData, '--', 'Color', groupColors{g}, 'LineWidth', 2)
        xlabel('Quartile')
        ylabel('Accuracy')
        ylim([0, 1])
        legEntries = [legEntries, {[id ' group' num2str(g)]}, {'Model'}]; %#ok<AGROW>
    end
    legend([{''; ''}; legEntries(:)], 'Interpreter', 'none', 'Location', 'best')
    title(['Sess1 ' model_name], 'Interpreter', 'none')
    saveas(fig, fullfile(plotDir, [id '_validate_group.svg']))
    close(fig)

    % Group p(engaged) distributions (Supplementary Fig 10 panels)
    num_bins = 10;
    bin_edges = linspace(0, 1, num_bins + 1);
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
    fig = figure('Position', [100, 100, 400 * length(groupIds), 400], 'Visible', 'off');
    for gi = 1:length(groupIds)
        g = groupIds(gi);
        dists = zeros(sum(groups == g), num_bins);
        idx = find(groups == g);
        for j = 1:numel(idx)
            dists(j, :) = histcounts(allLatent{idx(j)}, bin_edges, 'Normalization', 'probability');
        end
        subplot(1, length(groupIds), gi);
        bar(bin_centers, mean(dists, 1));
        hold on;
        er = errorbar(bin_centers, mean(dists, 1), std(dists, 1) / sqrt(size(dists, 1)));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        hold off;
        xlabel('p(engaged)');
        ylabel('Average Normalized Frequency');
        switch g
            case 1, title('Het');
            case 2, title('WT');
            case 3, title('KO');
        end
        xlim([0 1]);
        xticks(bin_edges);
        xtickangle(45);
        ylim([0, 1]);
    end
    sgtitle([id ' day 1']);
    saveas(fig, fullfile(plotDir, [id '_p_engaged_day1.svg']))
    close(fig)
end
end

function r = simulate_once(model_name, theta, this_data, latent, pReward)
%SIMULATE_ONCE  One simulated session; returns the simulated rewards.
%   For the Shank3B deterministic cohorts the generative model is called
%   without a reward-probability argument, which makes reward deterministic
%   given correctness and consumes one fewer random draw per trial.
if isempty(pReward)
    sim_data = feval(model_name, theta, this_data, latent);
else
    sim_data = feval(model_name, theta, this_data, latent, pReward);
end
r = sim_data(:, 3);
end
