function R = compare_params(id, opts)
%COMPARE_PARAMS  Group comparison of winning-model (a0b1s_hybrid) parameters.
%
%   R = COMPARE_PARAMS(id) loads the saved model fits of cohort ID, extracts
%   the MAP parameters of the winning model (a0b1s_hybrid), and compares the
%   Het and WT groups with two-tailed Mann-Whitney (ranksum) tests. These
%   are the comparisons reported in Figs 2-4 and Supplementary Figs 2, 5
%   and 7 of the manuscript. For cohorts that include KO animals (Shank3B
%   deterministic), the KO animals are shown in the figures but are not
%   part of the Het-vs-WT test.
%
%   Transformations (as in the manuscript):
%     alpha+       -> log(alpha+)
%     bias         -> |bias - 0.5|   (normalized, un-sided bias magnitude)
%     alpha+ * beta-> log(alpha+ * beta) (post-hoc composite)
%
%   R = COMPARE_PARAMS(id, opts) with opts.makePlots = false skips figures.
%
%   The returned struct R contains the test results:
%     pnames       - model parameter names
%     pvals        - Het-vs-WT ranksum p-value per parameter
%     alphaBetaPval - p-value for log(alpha+ * beta)
%     groups       - group index per subject (1 = Het, 2 = WT, 3 = KO)
%     corrMatrix, corrPvals - Kendall correlation between parameters
%
%   NOTE (scope):
%     1. Group comparisons of the p(engaged) mean and variance are
%        intentionally omitted: the values reported in the manuscript
%        cannot be reproduced from the archived fits (they were generated
%        by analysis code outside this folder). The engagement analyses
%        provided here are the KL-divergence distribution tests in
%        P_ENGAGED_TESTS.
%     2. The bias comparison is omitted for the Shank3B female cohort: the
%        archived fits yield p = 1.0 (verified with MATLAB and scipy),
%        while the manuscript reports p = 0.93, which likewise appears to
%        come from analysis code outside this folder.
%
%   See also COHORTS, FIT_MODELS, VALIDATE_MODELS, P_ENGAGED_TESTS.

%   Author: Jing-Jing Li (jl3676@berkeley.edu)

arguments
    id (1, :) char
    opts.makePlots (1, 1) logical = true
end

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'models', 'likelihood'));

C = cohorts(id);
model_name = 'a0b1s_hybrid';

m = load(fullfile(repoRoot, 'results', C.fitFile));
X = m.X; subjects = m.subjects; Ms = m.Ms; All_Params = m.All_Params;
model_ind = strcmp(cellfun(@(x) x.name, Ms(:), 'UniformOutput', false), model_name);
model_params = All_Params{model_ind};

%% Assign animal groups
[groups, groupLabels] = assign_groups(X, subjects, C.groupScheme);
het = find(groups == 1);
wt = find(groups == 2);

%% Parameter transforms and group tests
pnames = Ms{model_ind}.pnames;
nParams = size(model_params, 2);
pvals = zeros(1, nParams);
pvalsExact = zeros(1, nParams);
hetVals = cell(1, nParams);
wtVals = cell(1, nParams);
for p = 1:nParams
    [hv, wv] = transform_param(model_params, pnames, p, het, wt);
    hetVals{p} = hv;
    wtVals{p} = wv;
    pvals(p) = ranksum(hv, wv);
    % The manuscript's p-values were finalized in GraphPad Prism, which
    % uses the exact Mann-Whitney method; also report exact p-values for
    % direct comparability (small 2nd-decimal differences from the
    % asymptotic values are expected, with identical conclusions).
    pvalsExact(p) = ranksum(hv, wv, 'method', 'exact');
end

% Post-hoc composite log(alpha+ * beta)
alpha_ind = strcmp(pnames, 'alpha+');
beta_ind = strcmp(pnames, 'beta');
hetAB = log(model_params(het, alpha_ind) .* model_params(het, beta_ind));
wtAB = log(model_params(wt, alpha_ind) .* model_params(wt, beta_ind));
alphaBetaPval = ranksum(hetAB, wtAB);
alphaBetaPvalExact = ranksum(hetAB, wtAB, 'method', 'exact');

% The bias comparison is omitted for Shank_F_det (see header note 3)
omitBias = strcmp(id, 'Shank_F_det');
if omitBias
    pvals(strcmp(pnames, 'bias')) = NaN;
end

% Pairwise Kendall correlations between parameters
[corrMatrix, corrPvals] = corr(model_params, 'Type', 'Kendall');

R = struct('cohort', id, 'pnames', {pnames}, 'pvals', pvals, 'pvalsExact', pvalsExact, ...
    'hetVals', {hetVals}, 'wtVals', {wtVals}, ...
    'alphaBetaPval', alphaBetaPval, 'alphaBetaPvalExact', alphaBetaPvalExact, ...
    'groups', groups, 'groupLabels', {groupLabels}, ...
    'corrMatrix', corrMatrix, 'corrPvals', corrPvals);

fprintf('[%s] Het-vs-WT (Mann-Whitney):\n', id);
for p = 1:nParams
    if isnan(pvals(p))
        fprintf('  %-8s p = n/a (omitted, see header note)\n', pnames{p});
    else
        fprintf('  %-8s p = %.4f\n', pnames{p}, pvals(p));
    end
end
fprintf('  %-8s p = %.4f\n', 'log(a*b)', alphaBetaPval);

if opts.makePlots
    make_figures(id, C, model_name, Ms{model_ind}, model_params, ...
        groups, groupLabels, pvals, alphaBetaPval, corrMatrix, corrPvals, repoRoot);
end
end

function [hv, wv] = transform_param(model_params, pnames, p, het, wt)
pname = pnames{p};
if contains(pname, 'alpha')
    hv = log(model_params(het, p));
    wv = log(model_params(wt, p));
elseif strcmp(pname, 'bias')
    hv = abs(model_params(het, p) - 0.5);
    wv = abs(model_params(wt, p) - 0.5);
else
    hv = model_params(het, p);
    wv = model_params(wt, p);
end
end

function [groups, labels] = assign_groups(X, subjects, scheme)
%ASSIGN_GROUPS  Group index per subject: 1 = Het, 2 = WT, 3 = KO.
switch scheme
    case 'two'
        labels = {'HT', 'WT'};
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
        labels = {'HT', 'WT', 'KO'};
        groups = zeros(length(subjects), 1);
        for k = 1:length(subjects)
            T = find((X.Animal_ID == subjects(k)) & X.schedule < 3);
            if X.genotype(T(1)) == 1
                groups(k) = 1;
            elseif X.genotype(T(1)) == 2
                groups(k) = 3;
            else
                groups(k) = 2;
            end
        end
end
end

function make_figures(id, C, model_name, fit_model, model_params, ...
    groups, groupLabels, pvals, alphaBetaPval, corrMatrix, corrPvals, repoRoot)

colors = {'#80B3FF', '#F0539F'};
groupIds = unique(groups)';

pnames = fit_model.pnames;
nParams = size(model_params, 2);
nPanels = nParams + 1; % + log(alpha+ * beta)
ncols = 5;
nrows = ceil(nPanels / ncols);

figure('Position', [300, 300, 1200, 200 * nrows], 'Visible', 'off');
for p = 1:nParams
    subplot(nrows, ncols, p)
    pname = pnames{p};
    if contains(pname, 'alpha')
        vals = log(model_params(:, p));
        ttl = ['log ' pname];
    elseif strcmp(pname, 'bias')
        vals = abs(model_params(:, p) - 0.5);
        ttl = '|bias - 0.5|';
    else
        vals = model_params(:, p);
        ttl = pname;
    end
    plot_group_bars(vals, groups, groupIds, groupLabels, colors)
    if ~contains(pname, 'alpha') && ~strcmp(pname, 'bias')
        ylim([fit_model.pMin(p) fit_model.pMax(p)])
    end
    if isnan(pvals(p))
        title(ttl)
    else
        title(sprintf('%s (p=%.3f)', ttl, pvals(p)))
    end
end

alpha_ind = strcmp(pnames, 'alpha+');
beta_ind = strcmp(pnames, 'beta');
subplot(nrows, ncols, nParams + 1)
plot_group_bars(log(model_params(:, alpha_ind) .* model_params(:, beta_ind)), groups, groupIds, groupLabels, colors)
title(sprintf('log beta*alpha (p=%.3f)', alphaBetaPval))

sgtitle([id ' ' model_name], 'Interpreter', 'none')

plotDir = fullfile(repoRoot, 'plots');
if ~exist(plotDir, 'dir'), mkdir(plotDir); end
saveas(gcf, fullfile(plotDir, [id '_compare_params.png']))
saveas(gcf, fullfile(plotDir, [id '_compare_params.svg']))
close(gcf)

% Pairwise parameter correlation matrix
figure('Visible', 'off');
imagesc(corrMatrix);
colorbar;
title(['Pairwise correlation between parameters (' model_name ')'], 'Interpreter', 'none');
xlabel('Columns'); ylabel('Columns');
axis square;
set(gca, 'XTick', 1:nParams, 'XTickLabel', pnames, 'YTick', 1:nParams, 'YTickLabel', pnames);
xtickangle(45);
for i = 1:size(corrMatrix, 1)
    for j = 1:size(corrMatrix, 2)
        if corrPvals(i, j) < 0.001
            sig_star = '***';
        elseif corrPvals(i, j) < 0.01
            sig_star = '**';
        elseif corrPvals(i, j) < 0.05
            sig_star = '*';
        else
            sig_star = '';
        end
        text(j, i, sprintf('%.2f%s', corrMatrix(i, j), sig_star), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
saveas(gcf, fullfile(plotDir, [id '_param_correlations.png']))
close(gcf)
end

function plot_group_bars(vals, groups, groupIds, labels, colors)
%Bar + swarm plot of one value per subject, grouped by genotype.
groupIds = groupIds(:)';
x = 1:length(groupIds);
ymean = arrayfun(@(g) median(vals(groups == g)), groupIds);
yerr = arrayfun(@(g) std(vals(groups == g)) / sqrt(sum(groups == g)), groupIds);
h = bar(x, ymean, 'FaceColor', colors{1}, 'EdgeColor', 'none');
hold on
er = errorbar(x, ymean, yerr, 'LineWidth', 1.5);
er.Color = [0 0 0];
er.LineStyle = 'none';
for gi = 1:length(groupIds)
    gv = vals(groups == groupIds(gi));
    swarmchart(repmat(h(1).XEndPoints(gi), numel(gv), 1), gv, 30, ...
        'MarkerFaceColor', colors{2}, 'MarkerEdgeColor', 'none', ...
        'XJitter', 'density', 'XJitterWidth', 0.3)
    hold on
end
xticks(x);
xticklabels(labels(groupIds));
end
