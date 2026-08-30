function outFile = fit_models(id, opts)
%FIT_MODELS  Fit all candidate models to Session-1 choice data of a cohort.
%
%   FIT_MODELS(id) fits the candidate models of cohort ID (see COHORTS) to
%   the cohort's Session-1 data with maximum a posteriori (MAP) estimation
%   and saves the results to results/<cohort fit file>. A normal prior
%   (mean 5, SD 7) is placed on the softmax inverse temperature beta; all
%   other parameters have uninformative priors. Optimization uses
%   GlobalSearch with an SQP fmincon local solver. For a0b3s_hybrid,
%   starting values of the shared parameters are initialized from the
%   a0b3s fits (a0b3s is therefore always fitted first).
%
%   FIT_MODELS(id, opts) allows optional fields:
%     subjectIndices - fit only these subject indices (default: all)
%     saveResults    - save the .mat results file (default: true)
%     resultsDir     - directory for the results file (default: repo results/)
%     makePlots      - make AIC/BIC model-comparison plots (default: true)
%     startValueSeed - if nonempty, seed the RNG as
%                      rng(k*1000+m+97*startValueSeed) before drawing each
%                      subject-model's starting values, making fits
%                      reproducible (by default, starting values are drawn
%                      from the unseeded RNG)
%
%   The saved file contains:
%     X          - behavior table (columns used: Animal_ID, genotype,
%                  schedule, action, reward1)
%     subjects   - subject IDs, in fitting order
%     Ms         - model specifications (see MODEL_SPECS)
%     All_Params - 1xM cell array; All_Params{m}(k, :) is the MAP parameter
%                  vector of subject k for model m
%     All_fits   - nSubjects x 6 x M array of fit measures:
%                  [subject index, neg log-likelihood, AIC, BIC, psr2, AIC0]
%     beta_mu, beta_sigma - prior hyperparameters
%
%   Requires: Optimization Toolbox, Global Optimization Toolbox,
%   Parallel Computing Toolbox, Statistics and Machine Learning Toolbox.
%
%   See also COHORTS, MODEL_SPECS, COMPARE_PARAMS, VALIDATE_MODELS.

%   Author: Jing-Jing Li (jl3676@berkeley.edu)

arguments
    id (1, :) char
    opts.subjectIndices = []
    opts.saveResults (1, 1) logical = true
    opts.resultsDir = ''
    opts.makePlots (1, 1) logical = true
    opts.startValueSeed = []
end

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'models', 'likelihood'));

C = cohorts(id);
d = load(fullfile(repoRoot, 'data', 'processed', C.dataFile));
X = d.X;
subjects = unique(X.Animal_ID);
if ~isempty(opts.subjectIndices)
    subjects = subjects(opts.subjectIndices);
end

beta_mu = 5;
beta_sigma = 7;

Ms = model_specs(C.models);

%% Fit models
All_Params = cell(length(Ms), 1);
All_fits = cell(length(Ms), 1);
for m = 1:length(Ms)
    fit_model = Ms{m};
    pmin = fit_model.pMin;
    pmax = fit_model.pMax;
    pdfs = fit_model.pdfs;

    fitmeasures = cell(length(subjects), 1);
    fitparams = cell(length(subjects), 1);

    parfor k = 1:length(subjects)
        s = subjects(k);
        this_data = get_subject_data(X, s);

        % Sample parameter starting values
        if ~isempty(opts.startValueSeed)
            rng(k * 1000 + m + 97 * opts.startValueSeed, 'twister'); % deterministic starting values (test mode)
        end
        par = zeros(length(pmin), 1);
        for p_ind = 1:length(pmin)
            par(p_ind) = pdfs{p_ind}(0); % 0 is an arbitrary placeholder
        end

        % Set starting values of hybrid model parameters to the best fit
        % parameters of the corresponding static model
        if contains(fit_model.name, 'a0b3s_hybrid')
            model_dynamic = fit_model;
            model_static_ind = find(strcmp(cellfun(@(x) x.name, Ms, 'UniformOutput', false), 'a0b3s'));
            model_static = Ms{model_static_ind};
            for z = 1:length(model_dynamic.pnames)
                this_p = model_dynamic.pnames{z};
                if sum(strcmp(this_p, model_static.pnames)) > 0
                    par(strcmp(this_p, model_dynamic.pnames)) = All_Params{model_static_ind}(k, strcmp(this_p, Ms{model_static_ind}.pnames));
                end
            end
        end

        % Objective function: negative log-likelihood + beta prior penalty
        llhfun = @(p) feval([fit_model.name, '_llh'], p, this_data);
        if sum(strcmp(fit_model.pnames, 'beta')) > 0
            beta_idx = strcmp(fit_model.pnames, 'beta');
            myfitfun = @(p) llhfun(p) + sum((p(beta_idx) - beta_mu).^2 ./ (2*beta_sigma.^2));
        else
            myfitfun = @(p) llhfun(p);
        end
        rng default % For reproducibility
        fmincon_opts = optimoptions(@fmincon, 'Algorithm', 'sqp');
        if C.inequalityConstraints
            % The probabilistic-cohort fits pass an empty
            % inequality-constraint set (numerically inert).
            Aineq = zeros(3, length(par));
            bineq = zeros(3, 1);
            problem = createOptimProblem('fmincon', 'objective', myfitfun, 'x0', par, 'lb', pmin, 'ub', pmax, ...
                'Aineq', Aineq, 'bineq', bineq, 'options', fmincon_opts);
        else
            problem = createOptimProblem('fmincon', 'objective', myfitfun, 'x0', par, 'lb', pmin, 'ub', pmax, 'options', fmincon_opts);
        end
        gs = GlobalSearch;
        [param, llh] = run(gs, problem);

        % Fit measures (AIC, BIC, pseudo-r2)
        ntrials = size(this_data, 1);
        AIC = 2 * llh + 2 * length(param);
        BIC = 2 * llh + log(ntrials) * length(param);
        AIC0 = -2 * log(1/3) * ntrials;
        psr2 = (AIC0 - AIC) / AIC0;

        fitmeasures{k} = [k llh AIC BIC psr2 AIC0];
        fitparams{k} = param';
    end

    All_Params{m} = cell2mat(fitparams);
    All_fits{m} = cell2mat(fitmeasures);
end

% Reformat All_fits matrix
temp = All_fits;
All_fits = zeros(length(subjects), size(temp{1}, 2), length(Ms));
for i = 1:length(Ms)
    All_fits(:, :, i) = temp{i};
end

if isempty(opts.resultsDir)
    opts.resultsDir = fullfile(repoRoot, 'results');
end
if opts.saveResults
    outFile = fullfile(opts.resultsDir, C.fitFile);
    save(outFile, 'X', 'subjects', 'Ms', 'All_Params', 'All_fits', 'beta_mu', 'beta_sigma');
else
    outFile = '';
end

if opts.makePlots
    plot_model_comparison(id, X, subjects, Ms, All_fits, repoRoot);
end
end

function plot_model_comparison(id, X, subjects, Ms, All_fits, repoRoot)
% Model-comparison plots: delta AIC and delta BIC relative to each
% subject's mean, plus signed-rank tests between the winning model
% (a0b1s_hybrid) and the a0b2s/a0b3s hybrid models (as reported in the
% manuscript, Fig 2B, Supplementary Figs 5A,D and 7A,D).
modelNames = cellfun(@(x) x.name, Ms, 'UniformOutput', false);
win = find(strcmp(modelNames, 'a0b1s_hybrid'));

for measure = {'AIC', 'BIC'}
    switch measure{1}
        case 'AIC', col = 3;
        case 'BIC', col = 4;
    end
    vals = squeeze(All_fits(:, col, :));
    mvals = vals - repmat(mean(vals, 2), 1, size(vals, 2));

    figure('Position', [300 300 900 400], 'Visible', 'off');
    subplot(1, 2, 1)
    hold on
    bar(mean(mvals))
    errorbar(mean(mvals), std(mvals) / sqrt(size(mvals, 1)))
    xticks(1:length(Ms));
    xticklabels(modelNames);
    set(gca, 'TickLabelInterpreter', 'none')
    ylabel(['\Delta ' measure{1}])
    title('Models')

    subplot(1, 2, 2)
    yline(0, '--')
    hold on
    if strcmp(measure{1}, 'BIC')
        % The two model comparisons reported in the manuscript
        other = find(strcmp(modelNames, 'a0b3s_hybrid'));
        [p, ~, z] = signrank(vals(:, other), vals(:, win));
        plot(sort(vals(:, other) - vals(:, win), 'descend'), '.', 'MarkerSize', 15)
        title('a0b1s vs. a0b3s hybrid')
        sgtitle(sprintf('%s (p=%.4f, z=%.3f)', id, p, z.zval));
    else
        other = find(strcmp(modelNames, 'a0b2s_hybrid'));
        p = signrank(vals(:, other), vals(:, win), 'tail', 'left');
        plot(sort(vals(:, other) - vals(:, win), 'descend'), '.', 'MarkerSize', 15)
        title('a0b1s vs. a0b2s hybrid')
        sgtitle(sprintf('%s (p=%.4f)', id, p));
    end
    ylabel(['\Delta ' measure{1}])
    xlabel('sorted participant')
    set(gca, 'fontsize', 14)

    plotDir = fullfile(repoRoot, 'plots');
    if ~exist(plotDir, 'dir'), mkdir(plotDir); end
    saveas(gcf, fullfile(plotDir, [id '_fit_sess1_' measure{1} '.png']))
    saveas(gcf, fullfile(plotDir, [id '_fit_sess1_' measure{1} '.svg']))
    close(gcf)
end
end
