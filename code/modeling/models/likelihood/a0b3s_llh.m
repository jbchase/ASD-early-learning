function llh = a0b3s_llh(theta,data)

beta = theta(1);
alpha = [1e-6 theta(2)];%.2;
stick = [theta(3); theta(4); theta(5); theta(5)];
epsilon = 1e-6;%theta(7);

stimuli = data(:,1);
choices = data(:,2);
rewards = data(:,3);

Q = [.5 .5;.5 .5];
lt = log(.5);
llh = lt;
s = stimuli(1);
b(2) = epsilon/2 + (1-epsilon)/(1+exp(beta*(Q(s,1)-Q(s,2))));
b(1) = 1-b(2);

for k = 2:length(choices)
    s = stimuli(k-1);
    choice = choices(k-1);
    r = rewards(k-1);
    
    Q(s,choice) = Q(s,choice) + alpha(r+1)*(r-Q(s,choice));
      
    side=zeros(4,1);
    if stimuli(k-1) == stimuli(k) && rewards(k-1) == 0
        side(1) = 2 * (1.5 - choices(k-1)); % 1 for A1 and -1 for A2
    end
    if stimuli(k-1) ~= stimuli(k) && rewards(k-1) == 0
        side(2) = 2 * (1.5 - choices(k-1));
    end
    if stimuli(k-1) == stimuli(k) && rewards(k-1) > 0
        side(3) = 2 * (1.5 - choices(k-1));
    end
    if stimuli(k-1) ~= stimuli(k) && rewards(k-1) > 0
        side(4) = 2 * (1.5 - choices(k-1));
    end

    b(2) = epsilon/2 + (1-epsilon)/(1+exp(beta*(Q(stimuli(k),1)-Q(stimuli(k),2)+sum(stick.*side))));
    b(1) = 1-b(2);
       
    lt = log(b(choices(k)));
    llh = llh + lt;
end

llh = -llh;

end