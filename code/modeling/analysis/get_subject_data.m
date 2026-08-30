function this_data = get_subject_data(X, s)
%GET_SUBJECT_DATA  Extract one subject's Session-1 trial matrix.
%
%   this_data = GET_SUBJECT_DATA(X, s) returns an nTrials x 3 matrix
%   [schedule, action, reward] for subject ID s, restricted to odor-learning
%   trials (schedule < 3) with valid actions (action > 0; this also removes
%   missed/aborted trials, which are coded as action = 0 or NaN).
%
%   Column conventions:
%     schedule - 1 = odor A, 2 = odor B
%     action   - 1 = left port, 2 = right port
%     reward   - true if the trial was rewarded

%   Author: Jing-Jing Li (jl3676@berkeley.edu)

T = find((X.Animal_ID == s) & X.schedule < 3);
this_data = [X.schedule(T) X.action(T) X.reward1(T) > 0];
this_data = this_data(this_data(:, 2) > 0, :);
end
