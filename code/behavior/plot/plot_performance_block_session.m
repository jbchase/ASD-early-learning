function plot_performance_block_session(data, strain, genotype, savefigpath, tlabel)

if strcmp(tlabel,'AB') | strcmp(tlabel, 'CD') | strcmp(tlabel, 'AB-CD')
    nPlot = 3;
elseif strcmp(tlabel, 'DC') | strcmp(tlabel, 'AB-DC')
    nPlot = 6;
elseif strcmp(tlabel, 'AB-3h')
    nPlot = 3;
end
if contains(tlabel, 'reversed') | contains(tlabel, 'notreversed')
    if contains(tlabel, 'AB') | contains(tlabel, 'CD')
        nPlot = 3;
    else
        nPlot= 6;
    end
end

% control for 1000 trials maximum
x_plot = 1:10;
data = data(:,x_plot,:);

figure;
sgtitle(['Performance in the ', tlabel, ' sessions in block-session'])

for pp = 1:nPlot
    subplot(2,3,pp)

    genotype_list = unique(genotype);
    if strcmp(strain, 'Cntnap2_KO')  % only look at WT and KO
        genotype_list = {'KO', 'WT'};
    elseif strcmp(strain, 'Shank3B')
        genotype_list={'HET', 'WT'};
    end
    color_list = { 'red', 'blue','magenta'};
    % % temporary color for JC revision
    % colors = struct;
    % colors.TSC2_det_male = [72, 161, 217] / 255;
    % colors.TSC2_det_female = [22, 111, 56] / 255;
    % colors.Shank3_det_male = [50, 133, 141] / 255;
    % colors.Shank3_det_female = [248, 151, 29] / 255;
    % colors.TSC2_prob_male =  [125, 83, 162] / 255;
    % colors.Shank3_prob_male = [237, 28, 36] / 255;
    % color_list = {colors.TSC2_prob_male, 'black'};
    % compare performance within each block
    hold on;

    nBlocks = size(data,2);

    % Preallocate p-values for each block
    pvals = nan(1, nBlocks);

    for blockIdx = 1:nBlocks
        % Get data for this block for each genotype
        data1 = data(strcmp(genotype, genotype_list{1}), blockIdx, pp);
        data2 = data(strcmp(genotype, genotype_list{2}), blockIdx, pp);

        % Remove NaNs
        data1 = data1(~isnan(data1));
        data2 = data2(~isnan(data2));

        % Only test if both groups have >1 value
        if length(data1)>1 && length(data2)>1
            pvals(blockIdx) = ranksum(data1, data2);  % Mann–Whitney U test
        end
    end

    % --- Multiple comparisons correction (Bonferroni) ---
    p_fdr = mafdr(pvals,'BHFDR',true);  % fdr


    for gg = 1:length(genotype_list)
        % plot with blocks that have more than 3 animals
        not_nan_counts = sum(~isnan(data(strcmp(genotype, genotype_list{gg}),:,pp)));

        % Find indices of columns with more than 3 NaNs
        columns_with_nans = find(not_nan_counts  >= 3);

        mean_perf = nanmean(data(strcmp(genotype, genotype_list{gg}),columns_with_nans,pp),1);
        if sum(strcmp(genotype, genotype_list{gg}))>1
            ste_perf = nanstd(data(strcmp(genotype, genotype_list{gg}),columns_with_nans,pp),1)/sum(strcmp(genotype, genotype_list{gg}));
        else
            ste_perf = nan(length(columns_with_nans),1);
        end
        %plot(mean_perf)
        hold on;
        errorbar(columns_with_nans, mean_perf, ste_perf, 'LineWidth', 2,'Color', color_list{gg})

        % plot each session in thinner lines
        tempData= data(strcmp(genotype, genotype_list{gg}),:,pp);
        for tt =1:size(tempData,1)
            plot(x_plot, tempData(tt,:), 'LineWidth', 0.5, 'Color', color_list{gg},'LineStyle',':','HandleVisibility', 'off')
        end

    end
    yl = ylim;
    for blockIdx = 1:nBlocks
        if p_fdr(blockIdx) < 0.05
            text(blockIdx, yl(2)-0.05, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 15)
        end
    end

    ylim([0, 1]);
    %xticks([1,2,3,4]);
    %legend(genotype_list)
    title(['Session ', num2str(pp)])
    ylabel('Performance');
    % mixed-linear model
    if ismember('HEM', genotype_list) % nlgn3
        Hetdata = data(strcmp(genotype, 'HEM'),columns_with_nans,pp);
        nHet = sum(~all(isnan(Hetdata), 2));
    elseif ismember('KO', genotype_list) & length(genotype_list)==3 % if it's cntnap, compare KO with WT
        Hetdata = data(strcmp(genotype, 'HET'),columns_with_nans,pp);
        nHet = sum(~all(isnan(Hetdata), 2));
        KOdata = data(strcmp(genotype, 'KO'), columns_with_nans, pp);
        nKO = sum(~all(isnan(KOdata), 2));
    elseif ismember('KO', genotype_list) & length(genotype_list)==2
        Hetdata = data(strcmp(genotype, 'KO'),columns_with_nans,pp);
        nHet = sum(~all(isnan(Hetdata), 2));
    else
        Hetdata = data(strcmp(genotype, 'HET'),columns_with_nans,pp);
        nHet = sum(~all(isnan(Hetdata), 2));
    end
    WTdata = data(strcmp(genotype, 'WT'),columns_with_nans,pp);
    nWT = sum(~all(isnan(WTdata), 2));

    % check if there is empty data (not enough sessions/animals)
    if ~isempty(WTdata) && ~isempty(Hetdata)
        %% FDA test (functional data anslysis, nonparametric permutation test)
        stats_FDA = fda_longitudinal_permtest(WTdata, Hetdata, 10000);
        if length(genotype_list) == 2
            % WTdata, Hetdata are [N_blocks x N_subjects]
            nBlock = size(WTdata, 2);
            nWTsub = size(WTdata, 1);
            nHETsub = size(Hetdata,1);

            % Flatten data (stack all values)
            dat = [WTdata(:); Hetdata(:)];

            % Group labels
            group = [repmat({'WT'}, numel(WTdata), 1);
                repmat({'HET'}, numel(Hetdata), 1)];

            % Block labels (repeat block index for each subject)
            block = [repelem((1:nBlock)', nWTsub);
                repelem((1:nBlock)', nHETsub)];

            % Subject labels (repeat subject ID for each block)
            subject = [repmat((1:nWTsub)', nBlock, 1);
                repmat((1:nHETsub)' + nWTsub, nBlock, 1)];
        elseif length(genotype_list) == 3
            dat = [WTdata(:); Hetdata(:);KOdata(:)];
            group = [repmat('WT')];
        end
        %% Perform mixed ANOVA
       %% ==============================
        % 2-Way Repeated Measures ANOVA
        % Between: Group
        % Within:  Block
        % ==============================
        
        % Ensure matrices are [Subjects x Blocks]
        nBlock = size(WTdata,2);
        
        % Combine groups
        Y = [WTdata; Hetdata];
        
        Group = [repmat({'WT'}, size(WTdata,1),1); ...
                 repmat({'HET'}, size(Hetdata,1),1)];
        
        Subject = (1:size(Y,1))';
        
        % Remove subjects with ANY missing block (required for RM-ANOVA)
        valid_subjects = all(~isnan(Y),2);
        Y = Y(valid_subjects,:);
        Group = Group(valid_subjects);
        Subject = Subject(valid_subjects);
        
        %% ---- Create wide table ----
        Twide = array2table(Y);
        
        for b = 1:nBlock
            Twide.Properties.VariableNames{b} = sprintf('B%d',b);
        end
        
        Twide.Group   = categorical(Group);
        Twide.Subject = categorical(Subject);
        
        %% ---- Within-subject design ----
        within = table((1:nBlock)', 'VariableNames',{'Block'});
        %within.Block = categorical(within.Block);
        
        %% ---- Fit repeated-measures model ----
        rm = fitrm(Twide, ...
            sprintf('B1-B%d ~ Group',nBlock), ...
            'WithinDesign',within);
        
        %% ---- Run ANOVA ----
        ranovatbl = ranova(rm,'WithinModel','Block');   % Block + Interaction
        between   = anova(rm);                          % Group main effect
        disp(ranovatbl);
        
        %% ---- Extract p-values (robust to MATLAB versions) ----
        rowNames = string(ranovatbl.Properties.RowNames);

        % Group main effect
        pGroup = ranovatbl.pValue(rowNames == "Group");
        pLearn = ranovatbl.pValue(rowNames == "(Intercept):Block");
        pInt = ranovatbl.pValue(rowNames == "Group:Block");



        % % --- Remove NaNs ---
        % validIdx = ~isnan(dat);
        % dat = dat(validIdx);
        % group = group(validIdx);
        % block = block(validIdx);
        % subject = subject(validIdx);
        % 
        % % --- Convert to table for fitlme ---
        % tblData = table(dat, categorical(group), block, categorical(subject), ...
        %     'VariableNames', {'Y', 'Group', 'Block', 'Subject'});
        % 
        % % Ensure WT is the reference genotype
        % tblData.Group = reordercats(tblData.Group, {'WT', 'HET'});
        % tblData.Block2 = tblData.Block.^2;
        % tblData.lBlock = log(tblData.Block);
        % % --- Fit linear mixed-effects model ---
        % % Block as random intercept (repeated measure within Subject)
        % %lme = fitlme(tblData, 'Y ~ Group*Block + Group * Block2 + (Block|Subject)');
        % 
        % % generalized linear model
        % tblData.NumCorrect = round(tblData.Y * 100);   % Convert proportion → count
        % TotalTrials = 100;  
        % tblData.Prob = tblData.NumCorrect / 100;
        % logisticFun = @(b,x) b(1) ./ (1 + exp(-b(2)*(x-b(3))));
        % glme = fitglme(tblData, ...
        %     'NumCorrect ~ Group*Block +  (1+Block|Subject)', ...
        %     'Distribution','Binomial','Link','logit','BinomialSize',100);



        % glme_2 = fitglme(tblData, ...
        %     'NumCorrect ~ Group*Block + Group*Block2 +  (1+Block+Block2|Subject)', ...
        %     'Distribution','Binomial','Link','logit','BinomialSize',100);
        % 
        % glme_log = fitglme(tblData, ...
        %     'NumCorrect ~  Group*lBlock + Group*Block + (lBlock|Subject)', ...
        %     'Distribution','Binomial','Link','logit','BinomialSize',100);

        % --- Display results ---
        % disp(anova(glme)); % Mixed ANOVA table
        % pTable = anova(glme);
        % pGroup = pTable.pValue(strcmp(pTable.Term, 'Group'));
        % pLearn = pTable.pValue(strcmp(pTable.Term, 'Block'));
        % pInt   = pTable.pValue(strcmp(pTable.Term, 'Group:Block'));
        
        % --- individual learning slopes
        %% Extract subject-specific slopes from numeric randomEffects
        % RE = randomEffects(glme);   % numeric vector
        % 
        % % Determine number of subjects
        % Nsubj = length(unique(tblData.Subject));
        % 
        % % Model has (1 + Block|Subject) → 2 random effects per subject (intercept, slope)
        % subjectSlopes = RE(2:2:end);  % every 2nd element = slope
        % 
        % %% Add fixed slope to get total slope
        % fixedSlope = fixedEffects(glme);
        % slopeBlock = fixedSlope(strcmp(glme.CoefficientNames,'Block'));  % fixed slope
        % totalSlopes = slopeBlock + subjectSlopes;   % fitted slope per subject
        
        %% Map subjects to genotype
        % subjectIDs = unique(tblData.Subject);
 % convert to cell array
        
        % Split slopes by genotype
        % WTslopes  = totalSlopes(strcmp(genotype,'WT'));
        % HETslopes = totalSlopes(strcmp(genotype,'HET'));
        % 
        % %% Compare WT vs HET slopes
        % [p,h,stats]= ranksum(WTslopes,HETslopes);
        % 
        % fprintf('Mean slope WT:  %.4f\n', mean(WTslopes));
        % fprintf('Mean slope HET: %.4f\n', mean(HETslopes));
        % fprintf('p-value: %.3g\n', p);
        % 
        % --- Show p-values on figure ---
        text(0.5, 0.25, ['P(geno): ', num2str(pGroup, '%.3g')], 'FontSize', 15)
        text(0.5, 0.18, ['P(learn): ', num2str(pLearn, '%.3g')], 'FontSize', 15)
        text(0.5, 0.11, ['P(genoxlearn): ', num2str(pInt, '%.3g')], 'FontSize', 15)
    end
    if pp==nPlot
        legend(genotype_list, 'Location', 'east','FontSize', 14)
        legend('Box','off')
        text(1, 1, ['WT ', num2str(nWT)], 'FontSize', 10 );
        if ismember('KO', genotype_list) & length(genotype_list)==2
            text(5, 1, ['KO ', num2str(nHet)], 'FontSize', 12 );
        else
            text(5, 1, ['HET ', num2str(nHet)], 'FontSize', 12 );
        end
    end
end



%% savefigs
print(gcf,'-dpng',fullfile(savefigpath,['Performance in the  ', tlabel, ' sessions in block-session']));    %png format
saveas(gcf, fullfile(savefigpath, ['Performance in the ', tlabel, ' sessions in block-session']), 'fig');
saveas(gcf, fullfile(savefigpath, ['Performance in the ', tlabel, ' sessions in block-session']),'svg');

% close all including anova windows
delete(findall(0, 'Type', 'figure'))

%% compute AUC to compare the performance difference
% normalized over the trial blocks, and chance level (0.5)

perf_AUC = squeeze(nanmean((data-0.5), 2));
figure
for gg = 1:length(genotype_list)
    % plot with blocks that have more than 3 animals

    mean_perf = nanmean(perf_AUC(strcmp(genotype, genotype_list{gg}),:),1);

    ste_perf = nanstd(perf_AUC(strcmp(genotype, genotype_list{gg}),:),1)/sum(strcmp(genotype, genotype_list{gg}));

    %plot(mean_perf)
    hold on;
    errorbar(1:nPlot, mean_perf, ste_perf, 'LineWidth', 2,'Color', color_list{gg})

    % plot each session in thinner lines
    tempData= perf_AUC(strcmp(genotype, genotype_list{gg}),:);
    for tt =1:size(tempData,1)
        plot(1:nPlot, tempData(tt,:), 'LineWidth', 0.5, 'Color', color_list{gg},'LineStyle',':','HandleVisibility', 'off')
    end

end
ylim([-0.5, 0.5]);
%xticks([1,2,3,4]);
%legend(genotype_list)
ylabel('Performance (AUC)');
xticks(0:nPlot);
xticklabels(0:nPlot);
xlim([0.75,nPlot+0.25]);
xlabel('Sessions')
title(['AUC for ', tlabel]);
legend(genotype_list, 'Location', 'southeast')
legend('Box','off')

% ANOVA
if ismember('HEM', genotype_list) % nlgn3
    Hetdata = perf_AUC(strcmp(genotype, 'HEM'),:)';
    nHet = sum(~all(isnan(Hetdata), 2));
elseif ismember('KO', genotype_list) && length(genotype_list)==3 % if it's cntnap, compare KO with WT
    Hetdata = perf_AUC(strcmp(genotype, 'HET'),:)';
    nHet = sum(~all(isnan(Hetdata), 2));
    KOdata = perf_AUC(strcmp(genotype, 'KO'), :)';
    nHet = sum(~all(isnan(KOdata), 2));
elseif ismember('KO', genotype_list) & length(genotype_list)==2
    Hetdata = perf_AUC(strcmp(genotype, 'KO'),:)';
    nHet = sum(~all(isnan(Hetdata), 2));
else
    Hetdata = perf_AUC(strcmp(genotype, 'HET'),:)';
    nHet = sum(~all(isnan(Hetdata), 2));
end
WTdata = perf_AUC(strcmp(genotype, 'WT'),:)';
nWT = sum(strcmp(genotype, 'WT'));

% check if there is empty data (not enough sessions/animals)
if ~isempty(WTdata) && ~isempty(Hetdata)
    if length(genotype_list) == 2
        dat = [WTdata(:); Hetdata(:)];

        % Create grouping variables
        group = [repmat({'WT'}, numel(WTdata), 1); repmat({'HET'}, numel(Hetdata), 1)]; % Group labels
        session = [repmat((1:nPlot)', size(WTdata, 2), 1); repmat((1:nPlot)', size(Hetdata, 2), 1)];   % Quantile labels
    elseif length(genotype_list) == 3
        dat = [WTdata(:); Hetdata(:);KOdata(:)];
        session = [repmat((1:nPlot)', size(WTdata, 2), 1);
            repmat((1:nPlot)', size(Hetdata, 2), 1);
            repmat((1:nPlot)', size(KOdata, 2), 1)];
    end
    % Perform 2-way ANOVA
    [p, tbl, stats] = anovan(dat, {group, session}, ...
        'model', 'interaction', ...
        'varnames', {'Group', 'Session'});
    % display p-value for group in the figure

    text(1.5, 0.5, ['P(geno): ', num2str(p(1), '%.3g')], 'FontSize', 20 )
    text(1.5, 0.45, ['P(reversal): ', num2str(p(2), '%.3g')], 'FontSize', 20 )
end

print(gcf,'-dpng',fullfile(savefigpath,['Performance AUC in the  ', tlabel, ' sessions in block-session']));    %png format
saveas(gcf, fullfile(savefigpath, ['Performance in the AUC', tlabel, ' sessions in block-session']), 'fig');
saveas(gcf, fullfile(savefigpath, ['Performance in the AUC', tlabel, ' sessions in block-session']),'svg');

% close all including anova windows
delete(findall(0, 'Type', 'figure'))



