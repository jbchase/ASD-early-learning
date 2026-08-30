function plot_nTrials(nAB_trials, nCD_trials, nDC_trials,genotype, savefigpath)

% plot number of trials performed under different odors

figure;
sgtitle(['Number of trials performed'])


genotype_list = unique(genotype);
subplot(2,2,1)
x_plot = [1:size(nAB_trials,2)];
for gg = 1:length(genotype_list)
    mean_trials = nanmean(nAB_trials(strcmp(genotype, genotype_list{gg}),:),1);
    if sum(strcmp(genotype, genotype_list{gg}))>1
        ste_trials = nanstd(nAB_trials(strcmp(genotype, genotype_list{gg}),:),1)/sum(strcmp(genotype, genotype_list{gg}));
    else
        ste_trials = zeros(1, length(x_plot));
    end
    %plot(mean_perf)
    hold on;
    errorbar(x_plot, mean_trials, ste_trials, 'LineWidth', 2)
end
title('AB trials performed')

subplot(2,2,2)
x_plot = [4:4+size(nCD_trials,2)-1];
for gg = 1:length(genotype_list)
    mean_trials = nanmean(nCD_trials(strcmp(genotype, genotype_list{gg}),:),1);
    if sum(strcmp(genotype, genotype_list{gg}))>1
        ste_trials = nanstd(nCD_trials(strcmp(genotype, genotype_list{gg}),:),1)/sum(strcmp(genotype, genotype_list{gg}));
    else
        ste_trials = zeros(1, length(x_plot));
    end
    %plot(mean_perf)
    hold on;
    errorbar(x_plot, mean_trials, ste_trials, 'LineWidth', 2)
end
title('CD trials performed')
xlim([1, size(nAB_trials,2)]);

subplot(2,2,3)
x_plot = [7:7+size(nDC_trials,2)-1];
for gg = 1:length(genotype_list)
    mean_trials = nanmean(nDC_trials(strcmp(genotype, genotype_list{gg}),:),1);
    if sum(strcmp(genotype, genotype_list{gg}))>1
        ste_trials = nanstd(nDC_trials(strcmp(genotype, genotype_list{gg}),:),1)/sum(strcmp(genotype, genotype_list{gg}));
    else
        ste_trials = zeros(1, length(x_plot));
    end
    %plot(mean_perf)
    hold on;
    errorbar(x_plot, mean_trials, ste_trials, 'LineWidth', 2)
end
title('DC trials performed')
xlim([1, size(nAB_trials,2)]);

%xticks([1,2,3,4]);
legend(genotype_list)
legend('box', 'off')
ylabel('Performance');


%% savefigs
print(gcf,'-dpng',fullfile(savefigpath,['Number of trials performed']));    %png format
saveas(gcf, fullfile(savefigpath, ['Number of trials performed']), 'fig');
%savefig(fullfile(savefigpath, 'Number of trials performed.fig'));
saveas(gcf, fullfile(savefigpath, ['Number of trials performed']),'svg');

close;
