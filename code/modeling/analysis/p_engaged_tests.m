function R = p_engaged_tests(opts)
%P_ENGAGED_TESTS  KL-divergence permutation tests on p(engaged) distributions.
%
%   R = P_ENGAGED_TESTS() runs the three group-level p(engaged)
%   distribution comparisons reported in Supplementary Fig 10 of the
%   manuscript, all on Session 1 of the deterministic task:
%     A) Tsc2 males:     Het vs WT  (KL = 0.129, p = 0.433)
%     B) Shank3B males:  Het vs WT  (KL = 0.264, p = 0.20)
%     C) Shank3B males:  KO vs WT   (KL = 0.241, p = 0.266)
%
%   Group-level distributions are built by binning each animal's trial-level
%   p(engaged) values into 10 bins and averaging across animals. The
%   observed KL divergence between the two group distributions is compared
%   against a permutation null distribution (10,000 permutations of animal
%   labels); the p-value is the fraction of permuted KL values >= observed.
%
%   P_ENGAGED_TESTS(opts) allows optional fields:
%     nPermutations - number of permutations (default: 10000)
%     seed          - rng seed for reproducibility (default: 0)
%     makePlots     - plot permutation null distributions (default: true)
%
%   R is a struct array with fields comparison, groups, kl, p.
%
%   Group labels are assigned in subject order, matching the analysis
%   described in the manuscript. Note that comparison (A) therefore gives
%   slightly different values than originally reported (KL = 0.129,
%   p = 0.433), because the original computation mis-assigned genotypes to
%   animals; the conclusion (no significant distribution difference) is
%   unchanged.
%
%   See also KL_DIVERGENCE_PERMUTATION_TEST, GET_BINNED_P_ENGAGED.

%   Author: Jing-Jing Li (jl3676@berkeley.edu)

arguments
    opts.nPermutations (1, 1) double = 10000
    opts.seed (1, 1) double = 0
    opts.makePlots (1, 1) logical = true
end

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'analysis'));
rng(opts.seed);

comparisons = {
    'TSC_M_det',   'Tsc2 males, Het vs WT',    [1 2]; ...
    'Shank_M_det', 'Shank3B males, Het vs WT', [1 2]; ...
    'Shank_M_det', 'Shank3B males, KO vs WT',  [3 2]; ...
    };

R = struct('comparison', {}, 'groups', {}, 'kl', {}, 'p', {});
for c = 1:size(comparisons, 1)
    id = comparisons{c, 1};
    name = comparisons{c, 2};
    groupPair = comparisons{c, 3};

    C = cohorts(id);
    pe = load(fullfile(repoRoot, 'results', C.pEngagedFile));   % allLatent
    m = load(fullfile(repoRoot, 'results', C.fitFile));         % X, subjects
    X = m.X; subjects = m.subjects;

    % Group index per subject: 1 = Het, 2 = WT, 3 = KO
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

    if opts.makePlots
        [kl_div, p_value, ~] = kl_divergence_permutation_test( ...
            pe.allLatent, groups, groupPair, opts.nPermutations, name);
        plotDir = fullfile(repoRoot, 'plots');
        if ~exist(plotDir, 'dir'), mkdir(plotDir); end
        saveas(gcf, fullfile(plotDir, ['kl_permutation_' regexprep(name, '[^A-Za-z0-9]+', '_') '.png']));
        close(gcf);
    else
        [kl_div, p_value, ~] = kl_divergence_permutation_test( ...
            pe.allLatent, groups, groupPair, opts.nPermutations, '');
    end
    fprintf('%s: observed KL divergence = %.3f, permutation p = %.4f\n', name, kl_div, p_value);

    R(c).comparison = name;
    R(c).groups = groupPair;
    R(c).kl = kl_div;
    R(c).p = p_value;
end
end
