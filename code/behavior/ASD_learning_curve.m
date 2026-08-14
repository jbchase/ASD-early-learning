%% Plotting learning curve with trials per odor

group=SHANK3B_M_WT_AB_Int2;

GROUP='HET';
Data_analyzed='CD Novel Odor Session 1';
graph_output=0; % 1 for yes, 0 for no
graph_index=1; % 1 for WT, 2 for HOMO, 3 for HET

if graph_index==1
    C='b';D='bo-';E='c';F='g';G='co--';
elseif graph_index==2
    C='k';D='ko-';E='y';F='c';G='yo--';
else
    C='r';D='ro-';E='m';F='y';G='mo--';
end


session_type=1; %1= AB sessions %2= AB refresher trials in ABXX sessions %3= novel odor learning in ABXX sessions
trial_type=2; %1= Each session from each animal MUST have exactly the same number of 2 odor presentations.
              %2= Same number of 2 odor presentation within each session & Different number of presentation/odor between sessions. 

RP_number1=300; %number of presentation per odor (AB or novel odor learning)
RP_number2=50; %number of presentaiton per odor in AB refresher trials in ABXX sessions.

% number of trials in each bin, sliding
sliding_bin1=20; sliding_bin2=50;

% number of trials in each bin, not sliding
bin1=30; bin2=50; bin3=100;

% number of blocks within a session
block_1=3; block_2=4; block_3=10;

% iteration of delay since last correct to be analyzed
delay_iteration=10;

% % number of total sessions per animal being analyzed
% sessionN=1;

%%
warning('off','MATLAB:unknownElementsNowStruc');
warning('off');

% mean_perf_trials_per_odor=nan(length(group),900);
% mean_odor_perf_bin1s=nan(length(group),900);
% mean_odor_perf_bin2s=nan(length(group),900);
% mean_odor_perf_bin1=nan(length(group),90);
% mean_odor_perf_bin2=nan(length(group),18);

    Perf_Cum=nan(length(group),2000);
    Perf_per_sliding_bin1=nan(length(group),2000);
    Perf_per_sliding_bin2=nan(length(group),2000);
    Perf_per_bin1=nan(length(group),100);
    Perf_per_bin2=nan(length(group),60);
    Perf_per_bin3=nan(length(group),40);
    Perf_per_block_1=nan(length(group),block_1);
    Perf_per_block_2=nan(length(group),block_2);
    Perf_per_block_3=nan(length(group),block_3);
    Perf_per_odor=nan(length(group),2000);
    Perf_per_odor_sliding_bin1=nan(length(group),2000);
    Perf_per_odor_sliding_bin2=nan(length(group),2000);
    Perf_per_odor_bin1=nan(length(group),100);
    Perf_per_odor_bin2=nan(length(group),60);


    repeat_1_trial=nan(length(group),2000);
    repeat_2_trial=nan(length(group),2000);
    repeat_3_trial=nan(length(group),2000);
    nrepeat_1_trial=nan(length(group),2000);
    nrepeat_2_trial=nan(length(group),2000);
    nrepeat_3_trial=nan(length(group),2000);

    first_n_iteration=10;


for i=1:length(group)
    i;
    eval(['filerun=' group{i}])

    Perf_per_odor_session=nan(length(filerun),2000);
    Perf_per_odor_sliding_bin1_session=nan(length(filerun),2000);
    Perf_per_odor_sliding_bin2_session=nan(length(filerun),2000);
    Perf_per_odor_bin1_session=nan(length(filerun),100);
    Perf_per_odor_bin2_session=nan(length(filerun),60);
    


    for j=1:1 %:length(filerun)
%         
        data=load(filerun{j});
        exper=data.exper;
%         j=j-1; %included or not depending on which sessions.
        counted_trial(i,j)=exper.odor_rlwm.param.countedtrial.value;
        full_schedule=exper.odor_rlwm.param.schedule.value(1:counted_trial(i,j));
        full_result=exper.odor_rlwm.param.result.value(1:counted_trial(i,j));
        full_odorpokedur=exper.odor_rlwm.param.odorpokedur.value(1:counted_trial(end));
        full_portside=exper.odor_rlwm.param.port_side.value(1:counted_trial(i,j));
        Perf(i,j)=exper.odor_rlwm.param.totalscore.value;
        
        if Perf<0.5
            filerun{j}
            error('Overall session performance below 0.5')
        end
        
        
        % exclude trials with no data (missing trial?)
        full_schedule_1=full_schedule(full_result~=0);
        full_portside_1=full_portside(full_result~=0);
        full_result_1=full_result(full_result~=0);

        % exclude miss trials and withdraw trials
        full_schedule_1=full_schedule_1(full_result_1~=1.0&full_result_1~=3);
        full_portside_1=full_portside_1(full_result_1~=1.0&full_result_1~=3);
        full_result_1=full_result_1(full_result_1~=1.0&full_result_1~=3);


        % select the trials to be analyzed depending on session type
        if session_type==1
           % for AB session
            unique_schedule=unique(full_schedule_1);
            
        elseif session_type==2
            % for AB refresher trials in ABXX sessions
        
            full_portside_1=full_portside_1(ismember(full_schedule_1,[1 2]));
            full_result_1=full_result(ismember(full_schedule,[1 2])); 
            full_schedule_1=full_schedule_1(ismember(full_schedule_1,[1 2]));
            unique_schedule=unique(full_schedule_1);
            
        else
            % for new odor learning entire session
            % check last 60AB trials performance >=0.7
            AB_result=floor(full_result_1(ismember(full_schedule_1,[1 2])));
            AB_last60perf=sum(AB_result(end-59:end)==1)/60;
%             
%             if AB_last60perf <0.7
%                filerun{j}
%                error('ab last 60 trial performance failed (<0.7)')
%                
%             end
                          
            %for new odor learning

            full_portside_1=full_portside_1(~ismember(full_schedule_1,[1 2]));
            full_result_1=full_result_1(~ismember(full_schedule_1,[1 2]));
            full_schedule_1=full_schedule_1(~ismember(full_schedule_1,[1 2]));
            unique_schedule=unique(full_schedule_1);
        end
        
        Novel_odor_counted_trial(i,j)=length(full_result_1);
        % cumulative overall session learning curve
        % trial history in rows, individual animal in columns
        for m=1:length(full_result_1)
            Perf_Cum(i,m)=sum(floor(full_result_1(1:m))==1)/m;
        end

        if graph_output==1
        figure(i+10*(graph_index-1));hold on;
        subplot(2,5,1); hold on;
        plot(Perf_Cum(i,1:length(full_result_1)),C)
        plot([1 length(full_result_1)],[0.5 0.5],':k')
        plot([1 length(full_result_1)],[0.7 0.7],':k')
        xlabel('trials'),xlim([0 length(full_result_1)])
        ylabel('% Correct'), ylim([0 1])
        title('Cumulative')
        
        figure(graph_index*100+1);hold on;
        subplot(2,6,i);hold on;
        plot(Perf_Cum(i,1:length(full_result_1)),C)
        plot([1 length(full_result_1)],[0.5 0.5],':k')
        plot([1 length(full_result_1)],[0.7 0.7],':k')
        xlabel('trials'),xlim([0 length(full_result_1)])
        ylabel('% Correct'), ylim([0 1])
        title(group{i}(1:6))
        end

        % learning curve per trial with sliding window 1
        for k=sliding_bin1:length(full_result_1)
            Perf_per_sliding_bin1(i,k-sliding_bin1+1)=sum(floor(full_result_1(k-sliding_bin1+1:k))==1)/sliding_bin1;
        end

        if graph_output==1
        figure(i+10*(graph_index-1));hold on;
        subplot(2,5,2); hold on;
        plot(Perf_per_sliding_bin1(i,1:length(full_result_1)-sliding_bin1+1),C)
        plot([1 length(full_result_1)],[0.5 0.5],':k')
        plot([1 length(full_result_1)],[0.7 0.7],':k')
        xlabel('trials'),xlim([0 length(full_result_1)])
        ylabel('% Correct'), ylim([0 1])
        title('sliding 20 trials')

        figure(graph_index*100+2);hold on;
        subplot(2,6,i); hold on;
        plot(Perf_per_sliding_bin1(i,1:length(full_result_1)-sliding_bin1+1),C)
        plot([1 length(full_result_1)],[0.5 0.5],':k')
        plot([1 length(full_result_1)],[0.7 0.7],':k')
        xlabel('trials'), xlim([0 length(full_result_1)])
        ylabel('% Correct'), ylim([0 1])
        title(group{i}(1:6))
        end

        % learning curve per trial with sliding window 2
        for k=sliding_bin2:length(full_result_1)
            Perf_per_sliding_bin2(i,k-sliding_bin2+1)=sum(floor(full_result_1(k-sliding_bin2+1:k))==1)/sliding_bin2;
        end

        if graph_output==1    
        figure(i+10*(graph_index-1));hold on;
        subplot(2,5,3);hold on;
        plot(Perf_per_sliding_bin2(i,1:length(full_result_1)-sliding_bin2+1),C)
        plot([1 length(full_result_1)],[0.5 0.5],':k')
        plot([1 length(full_result_1)],[0.7 0.7],':k')
        xlabel('trials'),xlim([0 length(full_result_1)])
        ylabel('% Correct'), ylim([0 1])
        title('sliding 50 trials')

        figure(graph_index*100+3);hold on;
        subplot(2,6,i);hold on;
        plot(Perf_per_sliding_bin2(i,1:length(full_result_1)-sliding_bin2+1),C)
        plot([1 length(full_result_1)],[0.5 0.5],':k')
        plot([1 length(full_result_1)],[0.7 0.7],':k')
        xlabel('trials'),xlim([0 length(full_result_1)])
        ylabel('% Correct'), ylim([0 1])
        title(group{i}(1:6))
        end

        % learning curve per trial with bin 1
        for k=1:floor(length(full_result_1)/bin1)
            Perf_per_bin1(i,k)=sum(floor(full_result_1(bin1*(k-1)+1:bin1*k))==1)/bin1;
        end

        if graph_output==1
        figure(i+10*(graph_index-1));hold on;
        subplot(2,5,4);hold on;
        plot(Perf_per_bin1(i,1:floor(length(full_result_1)/bin1)),C)
        plot([1 floor(length(full_result_1)/bin1)],[0.5 0.5],':k')
        plot([1 floor(length(full_result_1)/bin1)],[0.7 0.7],':k')
        xlabel('Bins'),xlim([0 ceil(length(full_result_1)/bin1)])
        ylabel('% Correct'), ylim([0 1])
        title('30 trials/bin')

        figure(graph_index*100+4);hold on;
        subplot(2,6,i);hold on;
        plot(Perf_per_bin1(i,1:floor(length(full_result_1)/bin1)),C)
        plot([1 floor(length(full_result_1)/bin1)],[0.5 0.5],':k')
        plot([1 floor(length(full_result_1)/bin1)],[0.7 0.7],':k')
        xlabel('Bins'),xlim([0 ceil(length(full_result_1)/bin1)])
        ylabel('% Correct'), ylim([0 1])
        title(group{i}(1:6))
        end

         % learning curve per trial with bin 2
        for k=1:floor(length(full_result_1)/bin2)
            Perf_per_bin2(i,k)=sum(floor(full_result_1(bin2*(k-1)+1:bin2*k))==1)/bin2;
        end

        if graph_output==1
        figure(i+10*(graph_index-1));hold on;
        subplot(2,5,5);hold on;
        plot(Perf_per_bin2(i,1:floor(length(full_result_1)/bin2)),C)
        plot([1 floor(length(full_result_1)/bin2)], [0.5 0.5],':k')
        plot([1 floor(length(full_result_1)/bin2)],[0.7 0.7],':k')
        xlabel('Bins'),xlim([0 ceil(length(full_result_1)/bin2)])
        ylabel('% Correct'), ylim([0 1])
        title('50 trials/bin')

        figure(graph_index*100+5);hold on;
        subplot(2,6,i);hold on;
        plot(Perf_per_bin2(i,1:floor(length(full_result_1)/bin2)),C)
        plot([1 floor(length(full_result_1)/bin2)], [0.5 0.5],':k')
        plot([1 floor(length(full_result_1)/bin2)],[0.7 0.7],':k')
        xlabel('Bins'), xlim([0 ceil(length(full_result_1)/bin2)])
        ylabel('% Correct'), ylim([0 1])
        title(group{i}(1:6))
        end

        % learning curve per trial with bin 3
        for k=1:floor(length(full_result_1)/bin3)
            Perf_per_bin3(i,k)=sum(floor(full_result_1(bin3*(k-1)+1:bin3*k))==1)/bin3;
        end

        if graph_output==1
        figure(i+10*(graph_index-1));hold on;
        subplot(2,5,6); hold on;
        plot(Perf_per_bin3(i,1:floor(length(full_result_1)/bin3)),D)
        plot([1 floor(length(full_result_1)/bin3)], [0.5 0.5],':k')
        plot([1 floor(length(full_result_1)/bin3)],[0.7 0.7],':k')
        xlabel('Bins'),xlim([0 ceil(length(full_result_1)/bin3)])
        ylabel('% Correct'), ylim([0 1])
        title('100 trials/bin')

        figure(graph_index*100+6);hold on;
        subplot(2,6,i);hold on;
        plot(Perf_per_bin3(i,1:floor(length(full_result_1)/bin3)),D)
        plot([1 floor(length(full_result_1)/bin3)], [0.5 0.5],':k')
        plot([1 floor(length(full_result_1)/bin3)],[0.7 0.7],':k')
        xlabel('Bins'), xlim([0 ceil(length(full_result_1)/bin3)])
        ylabel('% Correct'), ylim([0 1])
        title(group{i}(1:6))
        end


        % learning curve per trial with block 1
        for k=1:block_1
            Perf_per_block_1(i,k)=mean(floor(full_result_1(floor(length(full_result_1)/block_1)*(k-1)+1:floor(length(full_result_1)/block_1*k)))==1);
        end

        if graph_output==1
        figure(i+10*(graph_index-1));hold on;
        subplot(2,5,7);hold on;
        plot(Perf_per_block_1(i,1:block_1),D)
        plot([1 block_1],[0.5 0.5],':k')
        plot([1 block_1],[0.7 0.7],':k')
        xlabel('Blocks'),xlim([0 4])
        ylabel('% Correct'), ylim([0 1])
        title('3 blocks')

        figure(graph_index*100+7);hold on;
        subplot(2,6,i);hold on;
        plot(Perf_per_block_1(i,1:block_1),D)
        plot([1 block_1],[0.5 0.5],':k')
        plot([1 block_1],[0.7 0.7],':k')
        xlabel('Blocks'),xlim([0 4])
        ylabel('% Correct'), ylim([0 1])
        title(group{i}(1:6))
        end


         % learning curve per trial with block 2
        for k=1:block_2
            Perf_per_block_2(i,k)=mean(floor(full_result_1(floor(length(full_result_1)/block_2)*(k-1)+1:floor(length(full_result_1)/block_2*k)))==1);
        end

        if graph_output==1
        figure(i+10*(graph_index-1));hold on;
        subplot(2,5,8); hold on;
        plot(Perf_per_block_2(i,1:block_2),D)
        plot([1 block_2],[0.5 0.5],':k')
        plot([1 block_2],[0.7 0.7],':k')
        xlabel('Blocks'),xlim([0 5])
        ylabel('% Correct'), ylim([0 1])
        title('4 blocks')

        figure(graph_index*100+8); hold on;
        subplot(2,6,i);hold on;
        plot(Perf_per_block_2(i,1:block_2),D)
        plot([1 block_2],[0.5 0.5],':k')
        plot([1 block_2],[0.7 0.7],':k')
        xlabel('Blocks'),xlim([0 5])
        ylabel('% Correct'), ylim([0 1])
        title(group{i}(1:6))
        end

         % learning curve per trial with block 3
        for k=1:block_3
            Perf_per_block_3(i,k)=mean(floor(full_result_1(floor(length(full_result_1)/block_3)*(k-1)+1:floor(length(full_result_1)/block_3*k)))==1);
        end

        if graph_output==1
        figure(i+10*(graph_index-1));hold on;
        subplot(2,5,9); hold on
        plot(Perf_per_block_3(i,1:block_3),D)
        plot([1 block_3],[0.5 0.5],':k')
        plot([1 block_3],[0.7 0.7],':k')
        xlabel('Blocks'),xlim([0 11])
        ylabel('% Correct'), ylim([0 1])
        title('10 blocks')

        figure(graph_index*100+9);hold on;
        subplot(2,6,i);hold on;
        plot(Perf_per_block_3(i,1:block_3),D)
        plot([1 block_3],[0.5 0.5],':k')
        plot([1 block_3],[0.7 0.7],':k')
        xlabel('Blocks'),xlim([0 11])
        ylabel('% Correct'), ylim([0 1])
        title(group{i}(1:6))
        end

%         % performance for repeated trial (1 2 3 times)
%         for k=2:length(full_schedule_1)
%             if full_schedule_1(k)==full_schedule_1(k-1)
%                 repeat_1_trial(k)=full_result_1(k);
%                 if k>=3 & full_schedule_1(k-1)==full_schedule_1(k-2)
%                     repeat_2_trial(k)=full_result_1(k);
%                     if k>=4 & full_schedule_1(k-2)==full_schedule_1(k-3)
%                         repeat_3_trial(k)=full_result_1(k);
%                     end
%                 end
%             else
%                 nrepeat_1_trial(k)=full_result_1(k);
%                 if k>=3 & full_schedule_1(k-1) ~= full_schedule_1(k-2)
%                     nrepeat_2_trial(k)=full_result_1(k);
%                     if k>=4 & full_schedule_1(k-2) ~= full_schedule_1(k-3)
%                         nrepeat_3_trial(k)=full_result_1(k);
%                     end
%                 end
%             end
%         end
% 
%         repeat_1_trial=repeat_1_trial(~isnan(repeat_1_trial));
%         repeat_2_trial=repeat_2_trial(~isnan(repeat_2_trial));
%         repeat_3_trial=repeat_3_trial(~isnan(repeat_3_trial));
%         nrepeat_1_trial=nrepeat_1_trial(~isnan(nrepeat_1_trial));
%         nrepeat_2_trial=nrepeat_2_trial(~isnan(nrepeat_2_trial));
%         nrepeat_3_trial=nrepeat_3_trial(~isnan(nrepeat_3_trial));
%         
% 
%         % delay since last correct, ignore odor identity
%         correct_trial_idx=[];correct_trial_2rp_idx=[]; correct_trial_3rp_idx=[]; 
%         
% 
%         correct_trial_idx=find(floor(full_result_1)==1);
%         
%         %full session
%         for k=1:4
%             correct_trial_delay_idx=correct_trial_idx+k;
%             correct_trial_delay_idx=correct_trial_delay_idx(find(correct_trial_delay_idx<=length(full_result_1)));
%             delay_since_correct=full_result_1(correct_trial_delay_idx);
%             delay_since_correct_perf(i,j,k)=sum(floor(delay_since_correct)==1)/length(delay_since_correct);
%         end
% 
%        
%    
% 
%         %selected iteration for delay since last correct comparison
%         delay_first_1_since_correct=delay_1_since_correct(1:delay_iteration);
%         delay_first_2_since_correct=delay_2_since_correct(1:delay_iteration);
%         delay_first_3_since_correct=delay_3_since_correct(1:delay_iteration);
% 
%         delay_after_1_since_correct=delay_1_since_correct(delay_iteration+1:end);
%         delay_after_2_since_correct=delay_2_since_correct(delay_iteration+1:end);
%         delay_after_3_since_correct=delay_3_since_correct(delay_iteration+1:end);
% 
%         delay_first_1_since_correct_perf(i,j)=sum(floor(delay_first_1_since_correct)==1)/length(delay_first_1_since_correct);
%         delay_first_2_since_correct_perf(i,j)=sum(floor(delay_first_2_since_correct)==1)/length(delay_first_2_since_correct);
%         delay_first_3_since_correct_perf(i,j)=sum(floor(delay_first_3_since_correct)==1)/length(delay_first_3_since_correct);
% 
%         delay_after_1_since_correct_perf(i,j)=sum(floor(delay_after_1_since_correct)==1)/length(delay_after_1_since_correct);
%         delay_after_2_since_correct_perf(i,j)=sum(floor(delay_after_2_since_correct)==1)/length(delay_after_2_since_correct);
%         delay_after_3_since_correct_perf(i,j)=sum(floor(delay_after_3_since_correct)==1)/length(delay_after_3_since_correct);
% 
% 
%         % performanace after 2 and 3 correct trials, ignore odor identity
%         for k=2:length(correct_trial_idx)
%             correct_trial_2rp_idx(k)=correct_trial_idx(k)-correct_trial_idx(k-1);
%         end
%         correct_trial_2rp_idx(1)=nan; % to correct information
%         for k=2:length(correct_trial_2rp_idx)
%             correct_trial_3rp_idx(k)=correct_trial_2rp_idx(k)-correct_trial_2rp_idx(k-1);
%         end
%         correct_trial_3rp_idx(1)=nan;% to correct information
% 
%         delay_2_since_2_correct=full_result_1(correct_trial_idx(find(correct_trial_2rp_idx==1))+1<=length(full_result_1));
%         delay_3_since_3_correct=full_result_1(correct_trial_idx(find(correct_trial_3rp_idx==0))+1<=length(full_result_1));
%         
%         delay_2_since_2_correct_perf(i,j)=sum(floor(delay_2_since_2_correct)==1)/length(delay_2_since_2_correct);
%         delay_3_since_3_correct_perf(i,j)=sum(floor(delay_3_since_3_correct)==1)/length(delay_3_since_3_correct);
%     

        % codes to allow different same_odor_length from each session within
        % animal and between animals
        for k=1:2
            odor_length(k)=length(find(full_schedule_1==unique_schedule(k)));
        end
        
        if trial_type==1
            if session_type~=2
                same_odor_length=RP_number1;
            else
                same_odor_length=RP_number2;
            end
        else
            same_odor_length=min(odor_length);            
        end

        same_odor_result=nan(2,max(same_odor_length));

        % find the same_odor_perf, learning curve trials per odor
        for k=1:2
            same_odor_schedule_idx=find(full_schedule_1==unique_schedule(k));
            same_odor_result=floor(full_result_1(same_odor_schedule_idx))==1;
                
            % cumulative performanc per odor and sliding bins per odor
            for m=1:length(same_odor_result)
                same_odor_perf(k,m)=sum(same_odor_result(1:m))/m;
                if m>=sliding_bin1
                    same_odor_perf_sliding_bin1(k,m-sliding_bin1+1)=sum(same_odor_result(m-sliding_bin1+1:m))/sliding_bin1;
                end
                if m>=sliding_bin2
                    same_odor_perf_sliding_bin2(k,m-sliding_bin2+1)=sum(same_odor_result(m-sliding_bin2+1:m))/sliding_bin2;
                end
            end

            % non-sliding bins per odor
            for n=1:floor(same_odor_length/bin1)
                same_odor_perf_bin1(k,n)=sum(same_odor_result((n-1)*bin1+1:bin1*n))/bin1;
            end
            for n=1:floor(same_odor_length/bin2)
                same_odor_perf_bin2(k,n)=sum(same_odor_result((n-1)*bin2+1:bin2*n))/bin2;
            end
        end

        % stored per odor performance and average per odor performance
        Perf_per_odor_session(j,1:same_odor_length)=mean(same_odor_perf(:,1:same_odor_length));
        Perf_per_odor_sliding_bin1_session(j,1:same_odor_length-sliding_bin1)=mean(same_odor_perf_sliding_bin1(:,1:same_odor_length-sliding_bin1));
        Perf_per_odor_sliding_bin2_session(j,1:same_odor_length-sliding_bin2)=mean(same_odor_perf_sliding_bin2(:,1:same_odor_length-sliding_bin2));
        Perf_per_odor_bin1_session(j,1:floor(same_odor_length/bin1))=mean(same_odor_perf_bin1(:,1:floor(same_odor_length/bin1)));
        Perf_per_odor_bin2_session(j,1:floor(same_odor_length/bin2))=mean(same_odor_perf_bin2(:,1:floor(same_odor_length/bin2)));
        
%         % delay since last correct, same odor
%         for k=1:2
%             correct_odor_idx=[]; correct_odor_2rp_idx=[];correct_odor_3rp_idx=[];
%             delay_1_odor_since_correct=[]; delay_2_odor_since_correct=[]; delay_3_odor_since_correct=[];
%             delay_2_odor_since_2_correct=[]; delay_3_odor_since_3_correct=[];
% 
%             same_odor_schedule_idx=find(full_schedule_1==unique_schedule(k));
%             same_odor_result=floor(full_result_1(same_odor_schedule_idx))==1;
% 
%             correct_odor_idx=find(floor(same_odor_result(1:same_odor_length))==1);
%             delay_1_odor_since_correct=same_odor_result(find(correct_odor_idx+1<=same_odor_length));
%             delay_2_odor_since_correct=same_odor_result(find(correct_odor_idx+2<=same_odor_length));
%             delay_3_odor_since_correct=same_odor_result(find(correct_odor_idx+3<=same_odor_length));
% 
%             delay_1_odor_since_correct_perf_per_odor(k)=sum(delay_1_odor_since_correct)/length(delay_1_odor_since_correct);
%             delay_2_odor_since_correct_perf_per_odor(k)=sum(delay_2_odor_since_correct)/length(delay_2_odor_since_correct);
%             delay_3_odor_since_correct_perf_per_odor(k)=sum(delay_3_odor_since_correct)/length(delay_3_odor_since_correct);
% 
%             % performance after 2 and 3 correct trials, same odor
%             for m=2:length(correct_odor_idx)
%                 correct_odor_2rp_idx(m)=correct_odor_idx(m)-correct_odor_idx(m-1);
%             end
%             correct_odor_2rp_idx(1)=nan;
%             for m=2:length(correct_odor_2rp_idx)
%                 correct_odor_3rp_idx(m)=correct_odor_2rp_idx(m)-correct_odor_2rp_idx(m-1);
%             end
%             correct_odor_3rp_idx(1)=nan;
% 
%             delay_2_odor_since_2_correct=same_odor_result(correct_odor_idx(find(correct_odor_2rp_idx==1))+1<=same_odor_length);
%             delay_3_odor_since_3_correct=same_odor_result(correct_odor_idx(find(correct_odor_3rp_idx==0))+1<=same_odor_length);
% 
%             delay_2_odor_since_2_correct_perf_per_odor(k)=sum(delay_2_odor_since_2_correct)/length(delay_2_odor_since_2_correct);
%             delay_3_odor_since_3_correct_perf_per_odor(k)=sum(delay_3_odor_since_3_correct)/length(delay_3_odor_since_3_correct);
%         end
% 
%         % delay since last correct, same odor, per animal
%         delay_1_odor_since_correct_perf(i,j)=mean(delay_1_odor_since_correct_perf_per_odor);
%         delay_2_odor_since_correct_perf(i,j)=mean(delay_2_odor_since_correct_perf_per_odor);
%         delay_3_odor_since_correct_perf(i,j)=mean(delay_3_odor_since_correct_perf_per_odor);
%         delay_2_odor_since_2_correct_perf(i,j)=mean(delay_2_odor_since_2_correct_perf_per_odor);
%         delay_3_odor_since_3_correct_perf(i,j)=mean(delay_3_odor_since_3_correct_perf_per_odor);
% 
%             
        end
        
        %store per odor performance for each animal
        Perf_per_odor(i,:)=nanmean(Perf_per_odor_session);
        Perf_per_odor_sliding_bin1(i,:)=nanmean(Perf_per_odor_sliding_bin1_session);
        Perf_per_odor_sliding_bin2(i,:)=nanmean(Perf_per_odor_sliding_bin2_session);
        Perf_per_odor_bin1(i,:)=nanmean(Perf_per_odor_bin1_session);
        Perf_per_odor_bin2(i,:)=nanmean(Perf_per_odor_bin2_session);


        
        

    end
%    
   
%             % find mean performance trials per odor per mouse
%             mean_perf_trials_per_odor(i,:)=nanmean(odor_perf,1);
%             
%             % find mean performance trials per odor by bin (not sliding) per mouse
%             mean_odor_perf_bin1s(i,:)=nanmean(odor_perf_bin1s,1);
%             mean_odor_perf_bin2s(i,:)=nanmean(odor_perf_bin2s,1);
%             mean_odor_perf_bin1(i,:)=nanmean(odor_perf_bin1,1);
%             mean_odor_perf_bin2(i,:)=nanmean(odor_perf_bin2,1);
%             


% group average data and figures

Perf_trial_data={Perf_Cum,Perf_per_sliding_bin1, Perf_per_sliding_bin2, Perf_per_bin1, Perf_per_bin2,...
    Perf_per_bin3,Perf_per_block_1, Perf_per_block_2, Perf_per_block_3};

for n=1:length(Perf_trial_data)
    Perf_trial_data_group{n}=nanmean(Perf_trial_data{n});
    Perf_trial_data_group_sem{n}=nanstd(Perf_trial_data{n})./sqrt(sum(~isnan(Perf_trial_data{n})));
end

% group average graphs
if graph_output==1
for q=1:length(Perf_trial_data)
figure(1000+q);hold on
% if q>=3
for n=1:length(group)
    plot(Perf_trial_data{q}(n,:),E)
end
% end

if sum(isnan(Perf_trial_data_group{q}))~=0
    trialN_plotted=min(find(isnan(Perf_trial_data_group{q})))-1;
else
    trialN_plotted=length(Perf_trial_data_group{q});
end

Y1=Perf_trial_data_group{q}(1:trialN_plotted)+Perf_trial_data_group_sem{q}(1:trialN_plotted);
Y2=Perf_trial_data_group{q}(1:trialN_plotted)-Perf_trial_data_group_sem{q}(1:trialN_plotted);
h=fill([1:trialN_plotted fliplr(1:trialN_plotted)],[Y1(1:trialN_plotted) fliplr(Y2(1:trialN_plotted))],F);
set(h,'EdgeColor',[1 1 1],'FaceAlpha',0.4)
if sum(isnan(Perf_trial_data_group{q}))~=0
    plot(Perf_trial_data_group{q}(1:trialN_plotted),C,'linewidth',2)
else
    plot(Perf_trial_data_group{q}(1:trialN_plotted),D,'linewidth',2)
end

plot([1 trialN_plotted],[0.5 0.5],':k')
plot([1 trialN_plotted],[0.7 0.7],':k')
ylabel('% Correct'), ylim([0 1])
xlim([0 trialN_plotted+1])
end
end

%[0.4,0.4,0.4])%[0.6,0.6,0.6]) %[0.55, 0.55 0.55] light gray males
% set(h,'EdgeColor',[1 1 1],'FaceAlpha',0.2);


% % group average for repeat and non-repeat trial
% for q=1:length(group)
%     repeat_1_trial_perf(q)=(sum(repeat_1_trial(q,:)==1.2)+sum(repeat_1_trial(q,:)==1.3))/sum(~isnan(repeat_1_trial(q,:)));
%     repeat_2_trial_perf(q)=(sum(repeat_2_trial(q,:)==1.2)+sum(repeat_2_trial(q,:)==1.3))/sum(~isnan(repeat_2_trial(q,:)));
%     repeat_3_trial_perf(q)=(sum(repeat_3_trial(q,:)==1.2)+sum(repeat_3_trial(q,:)==1.3))/sum(~isnan(repeat_3_trial(q,:)));
%     nrepeat_1_trial_perf(q)=(sum(nrepeat_1_trial(q,:)==1.2)+sum(nrepeat_1_trial(q,:)==1.3))/sum(~isnan(nrepeat_1_trial(q,:)));
%     nrepeat_2_trial_perf(q)=(sum(nrepeat_2_trial(q,:)==1.2)+sum(nrepeat_2_trial(q,:)==1.3))/sum(~isnan(nrepeat_2_trial(q,:)));
%     nrepeat_3_trial_perf(q)=(sum(nrepeat_3_trial(q,:)==1.2)+sum(nrepeat_3_trial(q,:)==1.3))/sum(~isnan(nrepeat_3_trial(q,:)));
% end   
% 
% if graph_output==1
% figure(1100); hold on;
% subplot(1,2,1);hold on;
% for q=1:length(group)
%     plot([1 2 3],[repeat_1_trial_perf(q) repeat_2_trial_perf(q) repeat_3_trial_perf(q)],G)
% end
% errorbar([mean(repeat_1_trial_perf) mean(repeat_2_trial_perf) mean(repeat_3_trial_perf)], ...
%     [std(repeat_1_trial_perf)/sqrt(length(group)) std(repeat_2_trial_perf)/sqrt(length(group)) std(repeat_3_trial_perf)/sqrt(length(group))],C,'linewidth',1.5)
% ylabel('% Correct'), ylim([0 1])
% xlabel('repeated trial'),xlim([0 4])
% title('repeated trials')
% 
% subplot(1,2,2); hold on;
% for q=1:length(group)
%     plot([1 2 3],[nrepeat_1_trial_perf(q) nrepeat_2_trial_perf(q) nrepeat_3_trial_perf(q)],G)
% end
% errorbar([mean(nrepeat_1_trial_perf) mean(nrepeat_2_trial_perf) mean(nrepeat_3_trial_perf)], ...
%     [std(nrepeat_1_trial_perf)/sqrt(length(group)) std(nrepeat_2_trial_perf)/sqrt(length(group)) std(nrepeat_3_trial_perf)/sqrt(length(group))],C,'linewidth',1.5)
% ylabel('number of trials'), ylim([0 1])
% xlabel('not-repeated trial'), xlim([0 4])
% title('non repeated trials')
% 
% % group average for delay since last correct, either ignor odor identity of
% % same odor
% figure(1200); hold on;
% subplot(2,2,1);hold on;
% for q=1:length(group)
%     plot([1 2 3 4],[mean(delay_since_correct_perf(q,:,1)) mean(delay_since_correct_perf(q,:,2)) mean(delay_since_correct_perf(q,:,3)) ...
%         mean(delay_since_correct_perf(q,:,4))],G)
% end
% errorbar([mean(mean(delay_first_1_since_correct_perf)) mean(mean(delay_first_2_since_correct_perf)) mean(mean(delay_first_3_since_correct_perf))],...
%     [std(mean(delay_first_1_since_correct_perf,2))/sqrt(length(group)) std(mean(delay_first_2_since_correct_perf,2))/sqrt(length(group)) ...
%     std(mean(delay_first_3_since_correct_perf,2))/sqrt(length(group))],C,'linewidth', 1.5)
% xlabel('delay since last correct'), xlim([0 4]) 
% ylabel('% Correct'), ylim([0 1])
% title('no identity and history')
% 
% subplot(2,2,2);hold on;
% for q=1:length(group)
%     plot([1 2 3], [mean(delay_1_since_correct_perf(q,:)) mean(delay_2_since_2_correct_perf(q,:)) mean(delay_3_since_3_correct_perf(q,:))],G)
% end
% errorbar([mean(mean(delay_1_since_correct_perf,2)) mean(mean(delay_2_since_2_correct_perf)) mean(mean(delay_3_since_3_correct_perf))],...
%     [std(mean(delay_1_since_correct_perf,2))/sqrt(length(group)) std(mean(delay_2_since_2_correct_perf,2))/sqrt(length(group)) ...
%     std(mean(delay_3_since_3_correct_perf,2))/sqrt(length(group))],C,'linewidth',1.5)
% xlabel('delay since last correct'), xlim([0 4])
% ylabel('% Correct'), ylim([0 1])
% title('no identity but history')
% 
% subplot(2,2,3);hold on;
% for q=1:length(group)
%     plot([1 2 3], [mean(delay_1_odor_since_correct_perf(q,:)) mean(delay_2_odor_since_correct_perf(q,:)) mean(delay_3_odor_since_correct_perf(q,:))],G)
% end
% errorbar([mean(mean(delay_1_odor_since_correct_perf,2)) mean(mean(delay_2_odor_since_correct_perf,2)) mean(mean(delay_3_odor_since_correct_perf,2))],...
%     [std(mean(delay_1_odor_since_correct_perf,2))/sqrt(length(group)) std(mean(delay_2_odor_since_correct_perf,2))/sqrt(length(group))... 
%     std(mean(delay_3_odor_since_correct_perf,2))/sqrt(length(group))],C,'lineWidth',1.5)
% xlabel('delay since last correct'), xlim([0 4])
% ylabel('% Correct'), ylim([0 1])
% title('same odor but no history')
% 
% subplot(2,2,4); hold on;
% for q=1:length(group)
%     plot([1 2 3], [mean(delay_1_odor_since_correct_perf(q,:)) mean(delay_2_odor_since_2_correct_perf(q,:)) mean(delay_3_odor_since_3_correct_perf(q,:))], G)
% end
% errorbar([mean(mean(delay_1_odor_since_correct_perf,2)) mean(mean(delay_2_odor_since_2_correct_perf,2)) mean(mean(delay_3_odor_since_3_correct_perf,2))], ...
%     [std(mean(delay_1_odor_since_correct_perf,2))/sqrt(length(group)) std(mean(delay_2_odor_since_2_correct_perf,2))/sqrt(length(group)) ...
%     std(mean(delay_3_odor_since_3_correct_perf,2))/sqrt(length(group))],C,'linewidth',1.5)
% xlabel('delay since last correct'), xlim([0 4])
% ylabel('% Correct'), ylim([0 1])
% title('same odor and history')
% end

% give title for each figure
   if graph_output==1
       for q=1:length(group)
           figure((graph_index-1)*10+q);hold on;
           sgtitle([group{q}(1:6) ' ' GROUP ' ' Data_analyzed])
       end
       figure(i+graph_index);hold on;
       sgtitle([group{i}(1:6) GROUP ' ' Data_analyzed])
    
       figure(graph_index*100+1);hold on;
       sgtitle({['ASD TSC2 ' GROUP ' ' Data_analyzed],'Cumulative Learning Curve'})
    
       figure(graph_index*100+2);hold on;
       sgtitle({['ASD TSC2 ' GROUP ' ' Data_analyzed],'Learning Curve with sliding bin 20 trials'})
    
       figure(graph_index*100+3);hold on;
       sgtitle({['ASD TSC2 ' GROUP ' ' Data_analyzed],'Learning Curve with sliding bin 50 trials'})
    
       figure(graph_index*100+4);hold on;
       sgtitle({['ASD TSC2 ' GROUP ' ' Data_analyzed],'Learning Curve with 30 trials/bin'})
    
       figure(graph_index*100+5);hold on;
       sgtitle({['ASD TSC2 ' GROUP ' ' Data_analyzed],'Learning Curve with 50 trials/bin'})
    
       figure(graph_index*100+6);hold on;
       sgtitle({['ASD TSC2 ' GROUP ' ' Data_analyzed],'Learning Curve with 100 trials/bin'})
       
       figure(graph_index*100+7);hold on;
       sgtitle({['ASD TSC2 ' GROUP ' ' Data_analyzed],'Learning Curve with 3 blocks'})
    
       figure(graph_index*100+8);hold on;
       sgtitle({['ASD TSC2 ' GROUP ' ' Data_analyzed],'Learning Curve with 4 blocks'})
    
       figure(graph_index*100+9);hold on;
       sgtitle({['ASD TSC2 ' GROUP ' ' Data_analyzed],'Learning Curve with 10 blocks'})
   end

%%

   
   for q=1:9
       figure(1000+q);hold on;
       f=get(gca,'children');
       legend([f(32),f(18),f(3)],'WT','HOMO','HET')
   end


   figure(1001);hold on;
   xlabel('trials')
   title({['ASD SHANK3B ' Data_analyzed],'Cumulative Learning Curve'})

   figure(1002);hold on;
   sgtitle({['ASD SHANK3B ' Data_analyzed],'Learning Curve with sliding bin 20 trials'})

   figure(1003);hold on;
   sgtitle({['ASD SHANK3B ' Data_analyzed],'Learning Curve with sliding bin 50 trials'})

   figure(1004);hold on;
   sgtitle({['ASD SHANK3B ' Data_analyzed],'Learning Curve with 30 trials/bin'})

   figure(1005);hold on;
   sgtitle({['ASD SHANK3B ' Data_analyzed],'Learning Curve with 50 trials/bin'})

   figure(1006);hold on;
   sgtitle({['ASD SHANK3B ' Data_analyzed],'Learning Curve with 100 trials/bin'})

   figure(1007);hold on;
   sgtitle({['ASD SHANK3B ' Data_analyzed],'Learning Curve with 3 blocks'})

   figure(1008);hold on;
   sgtitle({['ASD SHANK3B ' Data_analyzed],'Learning Curve with 4 blocks'})

   figure(1009);hold on;
   sgtitle({['ASD SHANK3B ' Data_analyzed],'Learning Curve with 10 blocks'})

%    figure(1100);hold on;
%    for w=1:2
%        subplot(1,2,w); hold on;
%        plot([0 4],[0.5 0.5],':k')
%        plot([0 4],[0.7 0.7],':k')
%        f=get(gca,'Children');
%        legend([f(14),f(3)],'WT','HET')
%    end
%    sgtitle({['ASD TSC2 ' Data_analyzed],'repeated vs non repeated trial performance'})
% 
%    figure(1200);hold on;
%    for w=1:4
%        subplot(2,2,w); hold on;
%        plot([0 4],[0.5 0.5],':k')
%        plot([0 4],[0.7 0.7],':k')
%        f=get(gca,'Children');
%        legend([f(14),f(3)],'WT','HET')
%    end
%    sgtitle({['ASD TSC2 ' Data_analyzed],'Delay Since Last Correct performance'})
% 


%%
% trialN_plotted=400;
% 
% figure(200);hold on;
% h=fill([1:trialN_plotted fliplr(1:trialN_plotted)],[AB_HETp1(1:trialN_plotted) fliplr(AB_HETp2(1:trialN_plotted))],[0.4,0.4,0.4])%[0.6,0.6,0.6]) %[0.55, 0.55 0.55] light gray males
% set(h,'EdgeColor',[1 1 1],'FaceAlpha',0.2);
% plot(MEAN_TSC2_F_HETp_perf_trials_per_odor(1:trialN_plotted),'Color',[0.4,0.4,0.4],'linewidth',2)
% h=fill([1:trialN_plotted fliplr(1:trialN_plotted)],[AB_HETp1(1:trialN_plotted) fliplr(AB_HETp2(1:trialN_plotted))],'g')%[34,139,34]/255) dark green; [0,0.5,1] light blue
% set(h,'EdgeColor',[1 1 1],'FaceAlpha',0.2);
% plot(MEAN_TSC2_F_HETp_perf_trials_per_odor(1:trialN_plotted),'Color','g','linewidth',2)
% 
% xlabel('Trials per odor')
% ylabel('fraction correct')
% ylim([0 1])
% f=get(gca,'Children')
% legend([f(3),f(1)],'HET','HET')
% plot([0 trialN_plotted],[0.5 0.5],'k:')
% title('TSC2 F Probabilistic Early AB Learning')
% 
% set(gca,'FontSize',14,'FontName','Arial')

% %%
% figure(1); hold on;
% subplot(2,5,1);hold on;
% plot(per_trial_perf,'linewidth',2)
% xlabel('trials');ylabel('performance'),ylim([0 1])
% title('Learning curve (per trial)')
% 
% subplot(2,5,2);hold on;
% plot(per_trial_perf_bin1s,'linewidth',2)
% xlabel('trials');ylabel('performance'),ylim([0 1])
% title('Learning curve (sliding bin=10)')
% 
% subplot(2,5,3);hold on;
% plot(per_trial_per_bin2s,'linewidth',2)
% xlabel('trials');ylabel('performance'),ylim([0 1])
% title('Learning curve (sliding bin=50)')
% 
% subplot(2,5,4);hold on;
% plot(per_trial_perf_bin1,'linewidth',2)
% xlabel('trials');ylabel('performance'),ylim([0 1])
% title('Learning curve (bin avg=10)')
% 
% subplot(2,5,5);hold on;
% plot(per_trial_perf_bin2,'linewidth',2)
% xlabel('trials');ylabel('performance'),ylim([0 1])
% title('Learning curve (bin avg=50)')
% 
% subplot(2,5,6);hold on;
% plot(nanmean(odor_perf),'linewidth',2)
% xlabel('trials per odor');ylabel('performance'),ylim([0 1])
% title('Learning curve')
% 
% subplot(2,5,7);hold on;
% plot(nanmean(odor_perf_bin1s),'linewidth',2)
% xlabel('trials per odor');ylabel('performance'),ylim([0 1])
% title('Learning curve (sliding bin=10)')
% 
% subplot(2,5,8);hold on;
% plot(nanmean(odor_perf_bin2s),'linewidth',2)
% xlabel('trials per odor');ylabel('performance'),ylim([0 1])
% title('Learning curve (sliding bin=50)')
% 
% subplot(2,5,9);hold on;
% plot(nanmean(odor_perf_bin1),'linewidth',2)
% xlabel('trials per odor');ylabel('performance'),ylim([0 1])
% title('Learning curve (bin avg=10)')
% 
% subplot(2,5,10);hold on;
% plot(nanmean(odor_perf_bin2),'linewidth',2)
% xlabel('trials per odor');ylabel('performance'),ylim([0 1])
% title('Learning curve (bin avg=50)')