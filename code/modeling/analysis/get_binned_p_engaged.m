function binnedFrequencies = get_binned_p_engaged(p_engaged)
%GET_BINNED_P_ENGAGED  Bin trial-level p(engaged) values per animal.
%
%   binnedFrequencies = GET_BINNED_P_ENGAGED(p_engaged) takes a 1 x nAnimals
%   cell array of trial-level engagement probabilities and returns an
%   nAnimals x 10 matrix of histogram counts over 10 evenly spaced bins
%   covering [0, 1].
%
%   See also KL_DIVERGENCE_PERMUTATION_TEST, P_ENGAGED_TESTS.

%   Author: Jing-Jing Li (jl3676@berkeley.edu)

numBins = 10;
edges = linspace(0, 1, numBins + 1); % 11 edges for 10 bins from [0, 1]

binnedFrequencies = zeros(length(p_engaged), numBins);
for i = 1:length(p_engaged)
    binnedFrequencies(i, :) = histcounts(p_engaged{i}, edges);
end
end
