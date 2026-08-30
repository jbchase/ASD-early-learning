function [kl_div_original, p_value, permuted_kl_divs] = kl_divergence_permutation_test(p_engaged, group_labels, group_vals, n_permutations, name)
%KL_DIVERGENCE_PERMUTATION_TEST  Permutation test on group p(engaged) distributions.
%
%   [kl_div, p_value, permuted_kl_divs] = KL_DIVERGENCE_PERMUTATION_TEST(
%   p_engaged, group_labels, group_vals, n_permutations, name) computes the
%   Kullback-Leibler divergence between the group-averaged p(engaged)
%   distributions of the two groups indicated by GROUP_VALS and tests it
%   against a permutation null distribution.
%
%   Inputs:
%     p_engaged      - 1 x nAnimals cell array of trial-level p(engaged)
%                      trajectories
%     group_labels   - nAnimals x 1 group index per animal
%     group_vals     - the two group indices to compare, e.g. [1 2]
%     n_permutations - number of label permutations
%     name           - title for the summary figure (optional)
%
%   The group-level distribution is the across-animal average of binned
%   trial-level p(engaged) frequencies (10 bins over [0, 1], see
%   GET_BINNED_P_ENGAGED). The p-value is the fraction of permuted KL
%   divergences >= the observed divergence (one-tailed).
%
%   Each permutation is drawn fresh from the original labeling, matching
%   the permutation procedure described in the manuscript.
%
%   See also GET_BINNED_P_ENGAGED, P_ENGAGED_TESTS.

%   Author: Jing-Jing Li (jl3676@berkeley.edu)

kl_div_original = group_kl(p_engaged, group_labels, group_vals);

permuted_kl_divs = zeros(n_permutations, 1);
for perm = 1:n_permutations
    permuted_labels = group_labels(randperm(length(group_labels)));
    permuted_kl_divs(perm) = group_kl(p_engaged, permuted_labels, group_vals);
end

p_value = mean(permuted_kl_divs >= kl_div_original);

if nargin >= 5 && ~isempty(name)
    figure('Visible', 'off');
    histogram(permuted_kl_divs, 'Normalization', 'probability');
    hold on;
    xline(kl_div_original, 'r-', 'Observed KL Divergence');
    title(name);
    xlabel('KL Divergence');
    ylabel('Probability');
end
end

function kl_div = group_kl(p_engaged, group_labels, group_vals)
%GROUP_KL  KL divergence between the group-averaged binned distributions.
binned_1 = get_binned_p_engaged(p_engaged(group_labels == group_vals(1)));
binned_2 = get_binned_p_engaged(p_engaged(group_labels == group_vals(2)));

dist1 = mean(binned_1, 1);
dist2 = mean(binned_2, 1);

% Normalize to probability distributions
dist1 = dist1 / sum(dist1);
dist2 = dist2 / sum(dist2);

% Small epsilon to avoid log(0) or division by zero
epsilon = 1e-12;
dist1 = dist1 + epsilon;
dist2 = dist2 + epsilon;

% KL divergence: KL(group1 || group2)
kl_div = sum(dist1 .* log(dist1 ./ dist2));
end
