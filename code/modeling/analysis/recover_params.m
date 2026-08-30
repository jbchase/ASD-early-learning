function R = recover_params(id, opts)
%RECOVER_PARAMS  Parameter recovery analysis for the winning model.
%
%   R = RECOVER_PARAMS(id) runs the generate-and-recover analysis of cohort
%   ID for the winning model a0b1s_hybrid (manuscript Supplementary Fig 3):
%   choice data are simulated for each animal using its best-fit (MAP)
%   parameters as ground truth, given the model-estimated latent state
%   occupancy trajectory; the model is then fitted back to the simulated
%   data, and recovered parameters are compared to ground truth with
%   Spearman correlations.
%
%   RECOVER_PARAMS(id, opts) allows optional fields:
%     subjectIndices - recover only these subject indices (default: all)
%     fitFile        - load fits from this file instead of the cohort's
%                      archived results file
%     makePlots      - make the generate-recover figure (default: true)
%     startValueSeed - if true, seed the RNG as rng(it) before drawing each
%                      subject's starting values, making recovery
%                      reproducible
%
%   The returned struct R contains:
%     genrec     - 1 x nSubjects cell array of [true; recovered] parameters
%     latent_gen - 1 x nSubjects cell array of p(engaged) trajectories used
%                  to generate the simulated data
%     latent_rec - 1 x nSubjects cell array of p(engaged) trajectories from
%                  the recovered parameters
%     rho, pval  - Spearman correlation and p-value per parameter
%
%   Requires: Optimization Toolbox, Global Optimization Toolbox,
%   Parallel Computing Toolbox, Statistics and Machine Learning Toolbox.
%
%   See also COHORTS, FIT_MODELS.

%   Author: Jing-Jing Li (jl3676@berkeley.edu)

arguments
    id (1, :) char
    opts.subjectIndices = []
    opts.fitFile = ''
    opts.makePlots (1, 1) logical = true
    opts.startValueSeed (1, 1) logical = false
end

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'models', 'likelihood'));
addpath(fullfile(repoRoot, 'models', 'generative'));

C = cohorts(id);
model = 'a0b1s_hybrid';

if isempty(opts.fitFile)
    opts.fitFile = fullfile(repoRoot, 'results', C.fitFile);
end
m = load(opts.fitFile);
X = m.X; subjects = m.subjects; Ms = m.Ms; All_Params = m.All_Params;
% Rebuild the model specifications from source: the Ms structs archived in
% the results files contain function handles referencing the fitting
% machine's absolute paths, which are not portable. Definitions are
% identical (see MODEL_SPECS).
Ms = model_specs(cellfun(@(x) x.name, Ms(:), 'UniformOutput', false));
model_ind = find(strcmp(cellfun(@(x) x.name, Ms(:), 'UniformOutput', false), model));
fitted_params = All_Params{model_ind};
if ~isempty(opts.subjectIndices)
    subjects = subjects(opts.subjectIndices);
    fitted_params = fitted_params(opts.subjectIndices, :);
end

beta_mu = 5;
beta_sigma = 7;

names = Ms{model_ind}.pnames;
pmin = Ms{model_ind}.pMin;
pmax = Ms{model_ind}.pMax;
pdfs = Ms{model_ind}.pdfs;

n_iters = size(fitted_params, 1);
genrec = cell(n_iters, 1);
latent_gen = cell(n_iters, 1);
latent_rec = cell(n_iters, 1);

parfor it = 1:n_iters
    theta = fitted_params(it, :);

    % Initialize parameter values for optimization
    if opts.startValueSeed
        rng(it, 'twister'); % deterministic starting values (test mode)
    end
    par = zeros(length(pmin), 1);
    for p_ind = 1:length(pmin)
        par(p_ind) = pdfs{p_ind}(0);
    end

    this_data = get_subject_data(X, subjects(it));
    latent = feval([model, '_latent'], theta, this_data);
    latent_gen{it} = latent;

    % Simulate data based on the fitted parameters
    this_data = feval(model, theta, this_data, latent);

    % Fit the model to the simulated data (MAP, as in FIT_MODELS)
    llhfun = @(p) feval([model, '_llh'], p, this_data);
    if sum(strcmp(names, 'beta')) > 0
        beta_idx = strcmp(names, 'beta');
        myfitfun = @(p) llhfun(p) + sum((p(beta_idx) - beta_mu).^2 ./ (2*beta_sigma.^2));
    else
        myfitfun = @(p) llhfun(p);
    end
    rng default % For reproducibility
    fmincon_opts = optimoptions(@fmincon, 'Algorithm', 'sqp');
    problem = createOptimProblem('fmincon', 'objective', myfitfun, 'x0', par, 'lb', pmin, 'ub', pmax, 'options', fmincon_opts);
    gs = GlobalSearch;
    [param, llh] = run(gs, problem); %#ok<ASGLU>
    latent = feval([model, '_latent'], param, this_data);
    latent_rec{it} = latent;

    genrec{it} = [theta' param];
end

% Spearman correlation between ground-truth and recovered values
nParams = length(names);
rho = zeros(1, nParams);
pval = zeros(1, nParams);
for p = 1:nParams
    if contains(names{p}, 'alpha')
        x = cell2mat(cellfun(@(x) log(x(p, 1)), genrec, 'uni', 0));
        y = cell2mat(cellfun(@(x) log(x(p, 2)), genrec, 'uni', 0));
    else
        x = cell2mat(cellfun(@(x) x(p, 1), genrec, 'uni', 0));
        y = cell2mat(cellfun(@(x) x(p, 2), genrec, 'uni', 0));
    end
    [rho(p), pval(p)] = corr(x, y, 'Type', 'Spearman');
end

R = struct('cohort', id, 'model', model, 'pnames', {names}, ...
    'genrec', {genrec}, 'latent_gen', {latent_gen}, 'latent_rec', {latent_rec}, ...
    'rho', rho, 'pval', pval);

fprintf('[%s] parameter recovery (Spearman):\n', id);
for p = 1:nParams
    fprintf('  %-8s rho = %.3f, p = %.4f\n', names{p}, rho(p), pval(p));
end

if opts.makePlots
    nrows = 3;
    ncols = 3;
    figure('Position', [300 300 ncols*200 nrows*200], 'Visible', 'off');
    for p = 1:nParams
        subplot(nrows, ncols, p)
        if contains(names{p}, 'alpha')
            x = cell2mat(cellfun(@(x) log(x(p, 1)), genrec, 'uni', 0));
            y = cell2mat(cellfun(@(x) log(x(p, 2)), genrec, 'uni', 0));
        else
            x = cell2mat(cellfun(@(x) x(p, 1), genrec, 'uni', 0));
            y = cell2mat(cellfun(@(x) x(p, 2), genrec, 'uni', 0));
        end
        scatter(x, y, 30, '.')
        lsline
        hold on
        plot(xlim, xlim, 'r')
        title([names{p}, ' rho=', num2str(rho(p), 3), ', p=', num2str(pval(p), 3)], 'Interpreter', 'none')
        xlabel('True')
        ylabel('Recovered')
    end

    % log(alpha+ * beta) recovery
    if sum(strcmp(names, 'beta')) > 0
        alpha_ind = strcmp(names, 'alpha+');
        subplot(nrows, ncols, nParams + 1)
        x = cell2mat(cellfun(@(x) log(x(alpha_ind, 1) * x(strcmp(names, 'beta'), 1)), genrec, 'uni', 0));
        y = cell2mat(cellfun(@(x) log(x(alpha_ind, 2) * x(strcmp(names, 'beta'), 2)), genrec, 'uni', 0));
        scatter(x, y, 30, '.')
        lsline
        [rhoAB, pAB] = corr(x, y, 'Type', 'Spearman');
        title(['log(a*b) rho=', num2str(rhoAB, 3), ', p=', num2str(pAB, 3)], 'Interpreter', 'none')
        xlabel('True')
        ylabel('Recovered')
        hold on
        plot(xlim, xlim, 'r')
    end

    % p(engaged) mean and variance recovery
    latent_gen_mean = cellfun(@mean, latent_gen);
    latent_rec_mean = cellfun(@mean, latent_rec);
    subplot(nrows, ncols, nParams + 2)
    scatter(latent_gen_mean, latent_rec_mean, 30, '.')
    [rhoM, pM] = corr(latent_gen_mean', latent_rec_mean', 'Type', 'Spearman');
    lsline
    title(['p(e) mean rho=', num2str(rhoM, 3), ', p=', num2str(pM, 3)], 'Interpreter', 'none')
    xlabel('True')
    ylabel('Recovered')
    hold on
    plot(xlim, xlim, 'r')

    latent_gen_var = cellfun(@var, latent_gen);
    latent_rec_var = cellfun(@var, latent_rec);
    subplot(nrows, ncols, nParams + 3)
    scatter(latent_gen_var, latent_rec_var, 30, '.')
    [rhoV, pV] = corr(latent_gen_var', latent_rec_var', 'Type', 'Spearman');
    lsline
    title(['p(e) var rho=', num2str(rhoV, 3), ', p=', num2str(pV, 3)], 'Interpreter', 'none')
    xlabel('True')
    ylabel('Recovered')
    hold on
    plot(xlim, xlim, 'r')

    h = sgtitle([id ' ' model]);
    h.Interpreter = 'none';

    plotDir = fullfile(repoRoot, 'plots');
    if ~exist(plotDir, 'dir'), mkdir(plotDir); end
    saveas(gcf, fullfile(plotDir, [id '_genrec_' model '.png']))
    close(gcf)
end
end
