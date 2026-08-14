function engagement_analysis(data_session, cutoff, p_engaged_session, genotypes, genoMask, label, color_list, savefigpath)

data_sum = vertcat(data_session{:});
p_engaged_sum = vertcat(p_engaged_session{:});
% cutoff
if length(cutoff) == 1
    % cutoff threshold from pooled data
    data = data_sum(data_sum < cutoff);
    p_engaged = p_engaged_sum(data_sum < cutoff);
    genoMask = genoMask(data_sum < cutoff);
else
    cut_mask= cell(1, length(cutoff));
    for ss = 1:length(cutoff)
        cut_mask{ss} = data_session{ss} <= cutoff(ss);
    end
    cut_mask = vertcat(cut_mask{:});
    data = data_sum(cut_mask);
    p_engaged = p_engaged_sum(cut_mask);
    genoMask = genoMask(cut_mask);
    % cutoff threshold applied to sessions individually
end


engaged_HET = data(p_engaged>0.75 & genoMask==1);
disengaged_HET = data(p_engaged<0.25 & genoMask==1);
engaged_WT = data(p_engaged>0.75 & genoMask==0);
disengaged_WT = data(p_engaged<0.25 & genoMask==0);

% ANOVA
data_anova = [engaged_HET; disengaged_HET; engaged_WT; disengaged_WT];
group = [ones(size(engaged_HET));2*ones(size(disengaged_HET));
    3*ones(size(engaged_WT));4*ones(size(disengaged_WT)) ];
% Factor 1: Genotype (1 = HET, 2 = WT)
geno = [ones(size(engaged_HET));
        ones(size(disengaged_HET));
        2*ones(size(engaged_WT));
        2*ones(size(disengaged_WT))];

% Factor 2: Engagement (1 = engaged, 2 = disengaged)
engaged = [ones(size(engaged_HET));
           2*ones(size(disengaged_HET));
           ones(size(engaged_WT));
           2*ones(size(disengaged_WT))];

% save it to a csv file
% Convert numeric labels to strings
geno_label = strings(length(geno),1);
geno_label(geno==1) = "Mut";
geno_label(geno==2) = "WT";

engaged_label = strings(length(engaged),1);
engaged_label(engaged==1) = "engaged";
engaged_label(engaged==2) = "disengaged";

% Create table
T = table(data_anova, geno_label, engaged_label, ...
          'VariableNames', {'data','genotype','engagement'});

% Save to CSV
savecsvfile = fullfile(savefigpath(1:end-8), 'Results', [label,'ITL_geno_engaged.csv']);
writetable(T,savecsvfile);

figure;
title(label)
%subplot(1,2,1)

%data = [engaged_HET; disengaged_HET];
%group = [ones(size(engaged_HET));2*ones(size(disengaged_HET))];
violinplot(data_anova, group);
set(gca, 'YScale','log')
ylabel('Time (s)')
xticks([1 2 3 4])
xticklabels({'Engaged HET','Disengaged HET', 'Engaged WT', 'Disengaged WT'})
%title('HET')
[p1,~] = ranksum(engaged_HET, disengaged_HET);
[p2, ~] = ranksum(engaged_WT, disengaged_WT);
[p3,~] = ranksum(engaged_HET, engaged_WT);
[p4, ~] = ranksum(disengaged_HET, disengaged_WT);

pvals = [p1 p2 p3 p4];

% Step 2: multiple comparison correction
% Option a: Bonferroni
p_fdr = mafdr(pvals, 'BHFDR', true);

text(1.5, 20, sprintf('Mut(EvsD):p = %.3g', p_fdr(1)), ...
    'HorizontalAlignment', 'center', 'FontSize', 20);
text(3.5, 20, sprintf('WT(EvsD):p = %.3g', p_fdr(2)), ...
    'HorizontalAlignment', 'center', 'FontSize', 20);

text(2, 30, sprintf('E(WTvsMut):p = %.3g', p_fdr(3)), ...
    'HorizontalAlignment', 'center', 'FontSize', 20);
text(4, 30, sprintf('D(WTvsMut):p = %.3g', p_fdr(4)), ...
    'HorizontalAlignment', 'center', 'FontSize', 20);

savefitpath = fullfile(savefigpath,'latent');
print(gcf,'-dpng',fullfile(savefitpath, ['Engage_session_',label]));    %png format
saveas(gcf, fullfile(savefitpath, ['Engage_session_', label]), 'fig');
saveas(gcf, fullfile(savefitpath, ['Engage_session_', label]),'svg');

%% plot average rt in sessions 
% not doable - too few disengaged trials in some session

% genotypes: n×1 cell array ('WT' or 'mutGene')
% rt_average_engaged, rt_average_disengage: n×1 numeric vectors

isWT  = strcmp(genotypes, 'WT');
if sum(contains(genotypes, 'HET'))> 0
    mutGene = 'HET';
    isMut = strcmp(genotypes, 'HET');   % adjust name if needed
elseif sum(contains(genotypes, 'KO'))> 0
    mutGene = 'KO';
    isMut = strcmp(genotypes, 'KO'); 
end

nsess = length(data_session);
average_engaged = zeros(nsess,1);
average_disengaged = zeros(nsess,1);
for ii = 1:nsess
    average_engaged(ii) = nanmean(data_session{ii}(p_engaged_session{ii}>0.75 & data_session{ii}<cutoff));
    average_disengaged(ii) = nanmean(data_session{ii}(p_engaged_session{ii}<0.25 & data_session{ii}<cutoff));
end
positions = [1, 2, 4, 5]; % [WT_engage, WT_disengage, Mut_engage, Mut_disengage]

figure; hold on;
title(label)
% --- WT ---
b1 = boxchart(ones(sum(isWT),1)*positions(1), average_engaged(isWT), ...
    'BoxFaceColor', [0 0 0], 'BoxWidth', 0.6);
b2 = boxchart(ones(sum(isWT),1)*positions(2), average_disengaged(isWT), ...
    'BoxFaceColor', [0.6 0.6 0.6], 'BoxWidth', 0.6);

% --- Mutant ---
b3 = boxchart(ones(sum(isMut),1)*positions(3), average_engaged(isMut), ...
    'BoxFaceColor', [1 0.4 0.4], 'BoxWidth', 0.6);
b4 = boxchart(ones(sum(isMut),1)*positions(4), average_disengaged(isMut), ...
    'BoxFaceColor', [0.8 0 0], 'BoxWidth', 0.6);

% --- Axes formatting ---
xlim([0 6]);
ylabel('Average time');
set(gca, 'TickDir', 'out', 'Box', 'off');

% Label engaged/disengaged per genotype
xticks(positions);
xticklabels({'Engaged','Disengaged','Engaged','Disengaged'});
ax = gca;
ax.XTickLabelRotation = 15;

% Group labels ("WT" and "mutGene") below x-axis
text(mean(positions(1:2)), max(ylim)-0.05*range(ylim), 'WT', ...
    'HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold');
text(mean(positions(3:4)), max(ylim)-0.05*range(ylim), mutGene, ...
    'HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold');

% Optional dashed separator between genotypes
plot([3 3], ylim, 'k--', 'LineWidth', 1);

% --- Mann–Whitney tests ---
p1  = ranksum(average_engaged(isWT),  average_disengaged(isWT));
p2 = ranksum(average_engaged(isMut), average_disengaged(isMut));
p3 = ranksum(average_engaged(isWT), average_engaged(isMut));
p4 = ranksum(average_disengaged(isWT), average_disengaged(isMut));
pvals = [p1 p2 p3 p4];

% Step 2: multiple comparison correction
% Option a: Bonferroni

p_fdr = mafdr(pvals, 'BHFDR', true); 

yl = ylim;
y_offset = 0.08 * range(yl);

% WT p-value
%plot([positions(1) positions(2)], [1 1]*(yl(2)-y_offset), 'k', 'LineWidth', 1);
text(mean(positions(1:2)), yl(2)-2, sprintf('WT(DvsE): %.3f', p_fdr(1)), ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20);

% Mut p-value
%plot([positions(3) positions(4)], [1 1]*(yl(2)-y_offset), 'k', 'LineWidth', 1);
text(mean(positions(3:4)), yl(2)-2, sprintf('Mut(DvsE): %.3f', p_fdr(2)), ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20);

%plot([positions(1) positions(3)], [1 1]*(yl(2)-y_offset-0.05), 'k', 'LineWidth', 1);
text(mean([positions(1),positions(3)]), yl(2)-3, sprintf('E(WTvsMut): %.3f', p_fdr(3)), ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20);

%plot([positions(2) positions(4)], [1 1]*(yl(2)-y_offset-0.1), 'k', 'LineWidth', 1);
text(mean([positions(2),positions(4)]), yl(2)-4, sprintf('D(WTvsMut): %.3f', p_fdr(4)), ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20);

savefitpath = fullfile(savefigpath,'latent');
print(gcf,'-dpng',fullfile(savefitpath, ['engage_session', label]));    %png format
saveas(gcf, fullfile(savefitpath, ['engage_session', label]), 'fig');
saveas(gcf, fullfile(savefitpath, ['engage_session', label]),'svg');
