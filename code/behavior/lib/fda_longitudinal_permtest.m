function stats = fda_longitudinal_permtest(WTdata, HETdata, nperm)
% FDA permutation tests for learning curves
%
% INPUTS:
%   WTdata  : [nWT × nBlocks] performance (0–1)
%   HETdata : [nHET × nBlocks] performance (0–1)
%   nperm   : number of permutations (default = 10000)
%
% OUTPUT (stats struct):
%   p_group     : global offset difference
%   p_learn     : learning-dynamics (interaction) difference
%   p_transient : transient (Sup-T) difference
%   T_*         : observed statistics
%   muWT, muHET : group mean curves
%   t           : block index

if nargin < 3
    nperm = 10000;
end

WTdata  = double(WTdata);
HETdata = double(HETdata);

nWT  = size(WTdata,1);
nHET = size(HETdata,1);
nBlk = size(WTdata,2);

assert(size(HETdata,2) == nBlk, 'WT/HET block mismatch');

% -----------------------------
% Combine data
% -----------------------------
data  = [WTdata; HETdata];
group = [true(nWT,1); false(nHET,1)]; % true = WT

t = 1:nBlk;

% -----------------------------
% Smooth each subject (FDA step)
% -----------------------------
curves = nan(size(data));


for i = 1:size(data,1)
    yi = data(i,:);
    ok = ~isnan(yi);

    if sum(ok) < 3
        continue   % not enough points to smooth
    end

    sp = spaps(t(ok), yi(ok), 1e-4);

    % evaluate only at observed blocks
    curves(i, ok) = fnval(sp, t(ok));
end

% -----------------------------
% Observed group means
% -----------------------------
muWT  = nanmean(curves(group,:), 1);
muHET = nanmean(curves(~group,:), 1);

diff_mu = muWT - muHET;

% -----------------------------
% Observed statistics
% -----------------------------

% 1) Global offset (L2 distance)
T_group_obs = trapz(t, diff_mu.^2);

% 2) Learning dynamics (derivative difference)
dWT  = gradient(muWT);
dHET = gradient(muHET);
T_learn_obs = trapz(t, (dWT - dHET).^2);

% 3) Transient difference (Supremum)
T_trans_obs = max(abs(diff_mu));

% -----------------------------
% Permutation tests
% -----------------------------
T_group_perm     = zeros(nperm,1);
T_learn_perm     = zeros(nperm,1);
T_trans_perm     = zeros(nperm,1);

N = size(curves,1);

for p = 1:nperm
    gperm = group(randperm(N));

    mu1 = nanmean(curves(gperm,:), 1);
    mu2 = nanmean(curves(~gperm,:), 1);

    dmu = mu1 - mu2;

    % Offset
    T_group_perm(p) = trapz(t, dmu.^2);

    % Learning dynamics
    d1 = gradient(mu1);
    d2 = gradient(mu2);
    T_learn_perm(p) = trapz(t, (d1 - d2).^2);

    % Transient (Sup-T)
    T_trans_perm(p) = max(abs(dmu));
end

% -----------------------------
% P-values
% -----------------------------
stats = struct();

stats.p_group     = mean(T_group_perm >= T_group_obs);
stats.p_learn     = mean(T_learn_perm >= T_learn_obs);
stats.p_transient = mean(T_trans_perm >= T_trans_obs);

stats.T_group     = T_group_obs;
stats.T_learn     = T_learn_obs;
stats.T_transient = T_trans_obs;

stats.muWT = muWT;
stats.muHET = muHET;
stats.t = t;

end
