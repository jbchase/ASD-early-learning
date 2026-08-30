function newData = a0b1s_hybrid(theta,data,latent,prew)

beta = theta(1);
% beta = 8;
alpha = [1e-6 theta(2)];%.2;
% alpha = [1e-6 0.002];
stick = [0; 0; theta(3); theta(3)];
% stick = [0; 0; 0; 0];
epsilon = 1e-6;%theta(7);
lapse = theta(4);%theta(2);%.05;
recover = theta(5);%theta(3);%.1;
bias = theta(6);
% bias = 0.5;

% state 2: RL
% state 1: random
% lapse: probability of lapsing from 2 to 1
% ret: probability of returning from random to RL

stimuli = data(:,1);
choices = zeros(size(data,1), 1);
rewards = zeros(size(data,1), 1);
policy = zeros(size(data,1), 1);

Q = [.5 .5;.5 .5];
engaged = 1;

firstS = stimuli(1);

for t = 1:size(data,1)
    s = stimuli(t);

    side = zeros(4,1);
    if t > 1
        if stimuli(t-1) == s && rewards(t-1) == 0
            side(1) = 2 * (1.5 - choices(t-1)); % 1 for A1 and -1 for A2
        end
        if stimuli(t-1) ~= s && rewards(t-1) == 0
            side(2) = 2 * (1.5 - choices(t-1));
        end
        if stimuli(t-1) == s && rewards(t-1) > 0
            side(3) = 2 * (1.5 - choices(t-1));
        end
        if stimuli(t-1) ~= s && rewards(t-1) > 0
            side(4) = 2 * (1.5 - choices(t-1));
        end
    end

    % make a choice
    if nargin > 2
        engaged = rand < latent(t);
    end

%     engaged = 1;

    if engaged
        pr = 1 / (1 + exp(beta * (Q(s,1) - Q(s,2) + sum(stick.*side)))); % Probability of choosing action A2
        if rand < lapse
            engaged = 0; % Lapse with a probability
        end
    else
        mQ = mean(Q);
%         pr = 1 / (1 + exp(beta * (mQ(1) - mQ(2) + stick5))); % Probability of choosing action A2 during lapse state
        if firstS == 2
            pr = bias;
        else
            pr = 1 - bias;
        end
        if rand < recover
            engaged = 1; % Return to engaged state with a probability
        end
    end

    choice = 1 + (rand < pr); % Choose action A2 with probability pr
    if rand < epsilon
        choice = randsample([1 2], 1); % Exploration: Choose randomly between actions A1 and A2
    end

    correct = choice == s;
    r = correct;
    if nargin == 4
        r = (rand < prew) * correct;
    end

    if s == 2
        policy(t) = 1 / (1 + exp(beta * (Q(s,1) - Q(s,2) + sum(stick.*side))));
    else
        policy(t) = 1 - 1 / (1 + exp(beta * (Q(s,1) - Q(s,2) + sum(stick.*side))));
    end

    Q(s,choice) = Q(s,choice) + alpha(r+1) * (r - Q(s,choice));

    choices(t) = choice;
    rewards(t) = r;
end

newData = [stimuli choices rewards policy];
end