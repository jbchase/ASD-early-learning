function llh = RL_epsilon_llh(theta,data)

beta = theta(1);%.8;
lapse = 0;%theta(2);%.05;
ret = 1;%theta(3);%.1;
alpha = theta(2:3);%.2;
stick = theta(4);
epsilon = 1e-6; 
epsilon = epsilon + (1-epsilon)*theta(5);
% state 2: RL
% state 1: random
% lapse: probability of lapsing from 2 to 1
% ret: probability of returning from random to RL

T = [1-ret lapse;ret 1-lapse];

stimuli = data(:,1);
choices = data(:,2);
rewards = data(:,3);

Q = [.5 .5;.5 .5];
lt = log(.5);
llh = lt;
p = [0 1];
s = stimuli(1);
b(2) = epsilon/2 + (1-epsilon)/(1+exp(beta*(Q(s,1)-Q(s,2))));
b(1) = 1-b(2);


for k = 2:length(choices)
    s = stimuli(k-1);
    choice = choices(k-1);
    r = rewards(k-1);
    for h=0:1
        pnew(h+1) = 0.5*p(1)*T(h+1,1) + b(choice)*p(2)*T(h+1,2);
    end
    p = pnew/exp(lt);
    
    Q(s,choice) = Q(s,choice) + alpha(r+1)*(r-Q(s,choice));
      
    side=1;
    if choice==2
        side=-1;
    end
    b(2) = epsilon/2 + (1-epsilon)/(1+exp(beta*(Q(stimuli(k),1)-Q(stimuli(k),2)+stick*side)));
    b(1) = 1-b(2);
       
    lt = log(0.5*p(1) + b(choices(k))*p(2));
    llh = llh + lt;
end
% 
% vars = {'Att','CorSide','choice','correct','reward','Q(cor)'};
% figure
% for i=1:6
%     subplot(2,3,i)
% plot(data(:,i))
% title(vars{i})
% end
llh = -llh;

end