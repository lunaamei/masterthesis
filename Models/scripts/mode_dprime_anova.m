% =========================================================================
% UPDATED FOR CV!
% THIS IS THE LATEST VERSION!!!!
% RUN THE model_analysis_dprime_stats_v2.m SCRIPT FIRST
% Two-way ANOVA on d' selectivity values
% Factors:  condition (3 levels: natural, gatys, texform)  x  category (5 levels)
% Loop:     one ANOVA per ROI (3 ROIs)
% DV:       dprime_single_stimulus
%
% This script saves:
%
%   table_descriptives_condition_by_category.csv
%   table_anova_clean.csv
%   table_posthoc_main_condition.csv
%   table_posthoc_interaction.csv
%
%   all_rois_descriptives.csv
%   all_rois_anova.csv
%   all_rois_posthoc_main_condition.csv
%   all_rois_posthoc_interaction.csv
%
%   <model_name>_all_rois_descriptives_table.tex
%   <model_name>_all_rois_anova_table.tex
%
%   Run once for baseline model folder
%   Run once for DVD model folder
% =========================================================================

clear; clc; close all;

%% ===== USER SETTINGS =====
% Choose ONE model folder at a time:
baseOutDir = '/Users/lunameidoering/Desktop/CODE/Models/dvd_imagenet_cv/results/';
%baseOutDir = '/Users/lunameidoering/Desktop/CODE/Models/baseline_imagenet_cv/results';

% Set model name for export files
model_name = 'dvd';
%model_name = 'baseline';

alpha = 0.05;
yLims = [-1.5 3.5];

known_conditions = {'natural', 'gatys', 'texform'};
known_rois       = {'LOTC_body_left', 'LOTC_hand_left', 'LOTC_tool_left'};
known_categories = {'bodies','hands','tools','mani','nman'};

% ------------------------------------------------------------------
% Bonferroni correction factors — mirrors master_dprime_analysis_complete.m
%   nAnovaCorr:       same effect tested across 3 ROIs
%   nDescPosthocCorr: 3 pairwise condition comparisons per ROI x category cell
%   nBaseCorr:        3 ROIs x 5 categories per condition (one-sample vs zero)
%   nPairCorr:        3 conditions per ROI (highest vs second-highest)
% ------------------------------------------------------------------
nAnovaCorr       = 3;
nDescPosthocCorr = 3;
nBaseCorr        = 15;
nPairCorr        = 3;

%% ===== FIND AND PARSE ALL SUBFOLDERS =====
subfolders = dir(baseOutDir);
subfolders = subfolders([subfolders.isdir]);
subfolders = subfolders(~ismember({subfolders.name}, {'.','..','anova_results'}));

folder_info = struct('name', {}, 'condition', {}, 'roi', {});

for f = 1:numel(subfolders)
    fname = subfolders(f).name;

    matched_cond = '';
    for c = 1:numel(known_conditions)
        if contains(fname, known_conditions{c})
            matched_cond = known_conditions{c};
            break;
        end
    end

    matched_roi = '';
    for r = 1:numel(known_rois)
        if contains(fname, known_rois{r})
            matched_roi = known_rois{r};
            break;
        end
    end

    if isempty(matched_cond) || isempty(matched_roi)
        fprintf('Skipping (could not parse condition/ROI): %s\n', fname);
        continue;
    end

    folder_info(end+1).name    = fname; 
    folder_info(end).condition = matched_cond;
    folder_info(end).roi       = matched_roi;
end

fprintf('Parsed %d / %d folders successfully.\n', numel(folder_info), numel(subfolders));

%% ===== GLOBAL OUTPUT DIR =====
globalOutDir = fullfile(baseOutDir, 'anova_results');
if ~exist(globalOutDir, 'dir'); mkdir(globalOutDir); end

%% ===== COLLECTORS FOR COMBINED TABLES =====
all_desc_tables  = {};
all_anova_tables = {};
all_posthoc_cond = {};
all_posthoc_int  = {};
all_baseline_rows = {};   % one-sample t-tests vs zero (Bonferroni N = nBaseCorr = 15)
all_pair_rows     = {};   % highest vs second-highest category (Bonferroni N = nPairCorr = 3)

%% ===== LOOP OVER ROIs =====
for roi_idx = 1:numel(known_rois)

    current_roi = known_rois{roi_idx};

    fprintf('\n============================================================\n');
    fprintf('ROI: %s\n', current_roi);
    fprintf('============================================================\n');

    %% --- load data for this ROI ---
    roi_mask    = strcmp({folder_info.roi}, current_roi);
    roi_folders = folder_info(roi_mask);

    if numel(roi_folders) ~= numel(known_conditions)
        warning('Expected %d condition folders for ROI %s, found %d — skipping.', ...
            numel(known_conditions), current_roi, numel(roi_folders));
        continue;
    end

    all_tables = {};

    for f = 1:numel(roi_folders)
        csv_path = fullfile(baseOutDir, roi_folders(f).name, ...
            'dprime_selectivity_single_stimulus.csv');

        if ~isfile(csv_path)
            warning('Missing CSV: %s — skipping ROI.', csv_path);
            break;
        end

        T = readtable(csv_path);
        T.condition = repmat({roi_folders(f).condition}, height(T), 1);
        T.roi       = repmat({current_roi}, height(T), 1);
        all_tables{end+1} = T; 

        fprintf('  Loaded %d rows — condition: %s\n', height(T), roi_folders(f).condition);
    end

    if numel(all_tables) ~= numel(known_conditions)
        warning('Could not load all conditions for ROI %s — skipping.', current_roi);
        continue;
    end

    data = vertcat(all_tables{:});
    fprintf('  Total observations: %d\n', height(data));

    %% --- prepare variables ---
    dv        = data.dprime_single_stimulus;
    condition = categorical(data.condition, known_conditions);
    category  = categorical(data.category, known_categories);

    cats   = categories(category);
    n_cat  = numel(cats);
    n_cond = numel(known_conditions);

    %% --- per-ROI output folder ---
    outDir = fullfile(globalOutDir, current_roi);
    if ~exist(outDir, 'dir'); mkdir(outDir); end

    %% ============================================================
    %% DESCRIPTIVE STATISTICS
    %% ============================================================
    fprintf('\n  Cell means (condition x category):\n');
    fprintf('  %-15s  %-12s  %7s  %7s  %5s\n', 'condition', 'category', 'mean', 'SD', 'n');
    
    desc_rows = {};
    
    for j = 1:n_cat
        
        row = {current_roi, char(cats{j})};
        n_stimuli_ref = NaN;
    
        for i = 1:n_cond
            mask = condition == known_conditions{i} & category == cats{j};
            vals = dv(mask & ~isnan(dv));
    
            m  = round(mean(vals, 'omitnan'), 2);
            sd = round(std(vals,  'omitnan'), 2);
            n  = numel(vals);
    
            if i == 1
                n_stimuli_ref = n;
            end
    
            row = [row, {m, sd}]; 
        end
    
        % add n once
        row = [row, {n_stimuli_ref}];
    
        % --- condition post-hocs: raw + Bonferroni-corrected (N = nDescPosthocCorr = 3) ---
        v_nat = dv(category == cats{j} & condition == 'natural');
        v_gat = dv(category == cats{j} & condition == 'gatys');
        v_tex = dv(category == cats{j} & condition == 'texform');
    
        [~, p_ng] = ttest2(v_nat, v_gat);
        [~, p_nt] = ttest2(v_nat, v_tex);
        [~, p_gt] = ttest2(v_gat, v_tex);

        p_ng_bonf = min(p_ng * nDescPosthocCorr, 1);
        p_nt_bonf = min(p_nt * nDescPosthocCorr, 1);
        p_gt_bonf = min(p_gt * nDescPosthocCorr, 1);
    
        row = [row, {p_ng, p_ng_bonf, p_nt, p_nt_bonf, p_gt, p_gt_bonf}];
    
        desc_rows(end+1, :) = row; 
    end
    
    T_desc = cell2table(desc_rows, 'VariableNames', ...
        {'ROI','Category', ...
         'natural_mean','natural_sd', ...
         'gatys_mean','gatys_sd', ...
         'texform_mean','texform_sd', ...
         'n_stimuli', ...
         'p_NG','p_NG_bonf', ...
         'p_NT','p_NT_bonf', ...
         'p_GT','p_GT_bonf'});

    writetable(T_desc, ...
        fullfile(outDir, 'table_descriptives_condition_by_category.csv'), ...
        'Delimiter', ';');

    all_desc_tables{end+1} = T_desc; 

    %% ============================================================
    %% TWO-WAY ANOVA
    %% ============================================================
    fprintf('\n  === Two-way ANOVA ===\n');
    [p_anova, tbl, stats] = anovan( ...
        dv, ...
        {condition, category}, ...
        'model',    'interaction', ...
        'varnames', {'condition', 'category'}, ...
        'display',  'off');

    p_condition   = p_anova(1);
    p_category    = p_anova(2);
    p_interaction = p_anova(3);

   
    ss_condition   = cell2mat(tbl(2,2));
    ss_category    = cell2mat(tbl(3,2));
    ss_interaction = cell2mat(tbl(4,2));
    ss_error       = cell2mat(tbl(5,2));

    df_condition   = cell2mat(tbl(2,3));
    df_category    = cell2mat(tbl(3,3));
    df_interaction = cell2mat(tbl(4,3));
    df_error       = cell2mat(tbl(5,3));

    F_condition    = round(cell2mat(tbl(2,6)), 2);
    F_category     = round(cell2mat(tbl(3,6)), 2);
    F_interaction  = round(cell2mat(tbl(4,6)), 2);

    ss_total = ss_condition + ss_category + ss_interaction + ss_error;

    eta2_condition   = round(ss_condition   / ss_total, 3);
    eta2_category    = round(ss_category    / ss_total, 3);
    eta2_interaction = round(ss_interaction / ss_total, 3);

    % Raw p stored here; Bonferroni across nAnovaCorr ROIs added after combining
    anova_rows = {
        current_roi, 'condition',          F_condition,   df_condition,   df_error, p_condition,   eta2_condition;
        current_roi, 'category',           F_category,    df_category,    df_error, p_category,    eta2_category;
        current_roi, 'condition:category', F_interaction, df_interaction, df_error, p_interaction, eta2_interaction
        };

    T_anova_clean = cell2table(anova_rows, 'VariableNames', ...
        {'ROI','Effect','F','df1','df2','p','eta2'});

    writetable(T_anova_clean, fullfile(outDir, 'table_anova_clean.csv'), ...
        'Delimiter', ';');

    all_anova_tables{end+1} = T_anova_clean; 
    %% ============================================================
    %% POST-HOC 1: MAIN EFFECT OF CONDITION
    %% ============================================================
    cond_pairs = nchoosek(1:n_cond, 2);
    n_pairs    = size(cond_pairs, 1);
    cond_rows  = {};

    for r = 1:n_pairs
        c1 = known_conditions{cond_pairs(r,1)};
        c2 = known_conditions{cond_pairs(r,2)};
        v1 = dv(condition == c1);
        v2 = dv(condition == c2);

        [~, p_raw, ~, tst] = ttest2(v1, v2);
        p_bonf = min(p_raw * n_pairs, 1);

        cond_rows(end+1,:) = {current_roi, c1, c2, tst.tstat, p_raw, p_bonf}; 
    end

    T_cond = cell2table(cond_rows, 'VariableNames', ...
        {'ROI','Condition1','Condition2','t','p_raw','p_bonf'});

    writetable(T_cond, fullfile(outDir, 'table_posthoc_main_condition.csv'), ...
        'Delimiter', ';');

    all_posthoc_cond{end+1} = T_cond; 

    %% ============================================================
    %% POST-HOC 2: INTERACTION (CONDITION WITHIN EACH CATEGORY)
    %% ============================================================
    n_simple_pairs   = n_cat * n_pairs;
    interaction_rows = {};

    for j = 1:n_cat
        for r = 1:n_pairs
            c1 = known_conditions{cond_pairs(r,1)};
            c2 = known_conditions{cond_pairs(r,2)};
            v1 = dv(category == cats{j} & condition == c1);
            v2 = dv(category == cats{j} & condition == c2);

            [~, p_raw, ~, tst] = ttest2(v1, v2);
            p_bonf = min(p_raw * n_simple_pairs, 1);

            interaction_rows(end+1,:) = { ...
                current_roi, char(cats{j}), c1, c2, ...
                tst.tstat, p_raw, p_bonf}; 
        end
    end

    T_int = cell2table(interaction_rows, 'VariableNames', ...
        {'ROI','Category','Condition1','Condition2','t','p_raw','p_bonf'});

    writetable(T_int, fullfile(outDir, 'table_posthoc_interaction.csv'), ...
        'Delimiter', ';');

    all_posthoc_int{end+1} = T_int; 

    
    %% PLOTS: INTERACTION + CONDITION MAIN EFFECT
    %% Same style as fMRI script, but error bars = SEM across stimuli
    

    % Pretty category labels for legend
    pretty_categories = {'Bodies','Hands','Tools','Manip','Nman'};

    % --- compute interaction means and SEMs ---
    intMean = nan(n_cond, n_cat);
    intSEM  = nan(n_cond, n_cat);

    for g = 1:n_cat
        for i = 1:n_cond
            vals = dv(category == cats{g} & condition == known_conditions{i});
            vals = vals(~isnan(vals));

            if ~isempty(vals)
                intMean(i,g) = mean(vals, 'omitnan');
                intSEM(i,g)  = std(vals, 'omitnan') / sqrt(numel(vals));
            end
        end
    end

    % --- compute condition main-effect means and SEMs ---
    condMean = nan(1, n_cond);
    condSEM  = nan(1, n_cond);

    for i = 1:n_cond
        vals = dv(condition == known_conditions{i});
        vals = vals(~isnan(vals));

        if ~isempty(vals)
            condMean(i) = mean(vals, 'omitnan');
            condSEM(i)  = std(vals, 'omitnan') / sqrt(numel(vals));
        end
    end

    % Use the same color palette as the fMRI script
    catColors = [1, 0.5, 0; ...
                 0.3, 0.6, 0.3; ... 
                 0, 0, 0.5; ...
                 0.3, 0.75, 0.93; ...
                 1, 0, 1];

    % --- 1. Interaction Plot ---
    fig1 = figure('Position', [100 100 700 500], 'Visible', 'off');
    hold on;

    for g = 1:n_cat
        errorbar(1:n_cond, intMean(:,g), intSEM(:,g), '-o', ...
            'Color', catColors(g,:), ...
            'LineWidth', 2, ...
            'MarkerSize', 7, ...
            'MarkerFaceColor', catColors(g,:), ...
            'CapSize', 8);
    end

    ylim(yLims);
    

    set(gca, 'XTick', 1:n_cond, 'XTickLabel', known_conditions, 'FontSize', 14);
  
    ylabel("Average d'", 'FontSize', 14, 'FontWeight', 'bold');
    title(sprintf('%s: Interaction', strrep(current_roi, '_', ' ')), 'Interpreter', 'none');

    legend(pretty_categories, 'Location', 'bestoutside', 'Box', 'off');
    grid off; box off;

    exportgraphics(fig1, fullfile(outDir, sprintf('plot_%s_interaction.png', current_roi)), 'Resolution', 300);

    % --- 2. Condition Main Effect Plot ---
    fig2 = figure('Position', [100 100 600 450], 'Visible', 'off');
    hold on;

    errorbar(1:n_cond, condMean, condSEM, ...
        '-ko', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'CapSize', 10);

    
    ylim([-1 1]);

    set(gca, 'XTick', 1:n_cond, 'XTickLabel', known_conditions, 'FontSize', 11);
    xlabel('Image Condition', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel("Average d'", 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('%s: Condition Main Effect', strrep(current_roi, '_', ' ')), 'Interpreter', 'none');

    grid on; box off;

    exportgraphics(fig2, fullfile(outDir, sprintf('plot_%s_main_cond.png', current_roi)), 'Resolution', 300);

    close(fig1);
    close(fig2);

    %% ============================================================
    %% BASELINE T-TESTS (one-sample vs zero) AND PAIRWISE HIGHEST VS SECOND
    %% ============================================================
    for ci = 1:n_cond
        cond_name = known_conditions{ci};

        %% --- one-sample t-tests against zero ---
        for j = 1:n_cat
            vals = dv(category == cats{j} & condition == cond_name);
            vals = vals(~isnan(vals));

            m_val  = mean(vals, 'omitnan');
            sd_val = std(vals,  'omitnan');
            n_val  = numel(vals);

            if n_val >= 2
                [~, p_val, ~, tstats] = ttest(vals, 0);
                t_val  = tstats.tstat;
                df_val = tstats.df;
            else
                p_val  = NaN;  t_val = NaN;  df_val = NaN;
            end

            p_bonf_val = min(p_val * nBaseCorr, 1);

            all_baseline_rows(end+1, :) = { ...
                current_roi, cond_name, char(cats{j}), ...
                n_val, m_val, sd_val, t_val, df_val, p_val, p_bonf_val ...
            };
        end

        %% pairwise: highest vs second-highest category
        cell_means = nan(1, n_cat);
        for j = 1:n_cat
            vals = dv(category == cats{j} & condition == cond_name);
            cell_means(j) = mean(vals(~isnan(vals)), 'omitnan');
        end

        [~, ord] = sort(cell_means, 'descend', 'MissingPlacement', 'last');

        if numel(ord) >= 2 && ~any(isnan(cell_means(ord(1:2))))
            idx1 = ord(1);
            idx2 = ord(2);

            v1 = dv(category == cats{idx1} & condition == cond_name);
            v2 = dv(category == cats{idx2} & condition == cond_name);
            v1 = v1(~isnan(v1));
            v2 = v2(~isnan(v2));

            if numel(v1) >= 2 && numel(v2) >= 2
                [~, p_pair, ~, stats_pair] = ttest2(v1, v2);
                t_pair  = stats_pair.tstat;
                df_pair = round(stats_pair.df);   
                n_pair  = numel(v1);               % n of highest-category stimuli
            else
                p_pair = NaN;  t_pair = NaN;  df_pair = NaN;  n_pair = NaN;
            end

            p_pair_bonf = min(p_pair * nPairCorr, 1);

            all_pair_rows(end+1, :) = { ...
                current_roi, cond_name, ...
                char(cats{idx1}), char(cats{idx2}), ...
                n_pair, t_pair, df_pair, p_pair, p_pair_bonf ...
            }; 
        else
            all_pair_rows(end+1, :) = { ...
                current_roi, cond_name, 'n/a', 'n/a', ...
                NaN, NaN, NaN, NaN, NaN ...
            }; 
        end
    end


end



%% ============================================================
%% SAVE COMBINED CSV TABLES
%% ============================================================
if ~isempty(all_desc_tables)
    T_desc_all = vertcat(all_desc_tables{:});
    writetable(T_desc_all, fullfile(globalOutDir, 'all_rois_descriptives.csv'), ...
        'Delimiter', ';');
else
    T_desc_all = table();
end

if ~isempty(all_anova_tables)
    T_anova_all = vertcat(all_anova_tables{:});
    % Bonferroni correction across nAnovaCorr ROIs, applied per effect type
    T_anova_all.p_bonf = min(T_anova_all.p * nAnovaCorr, 1);
    writetable(T_anova_all, fullfile(globalOutDir, 'all_rois_anova.csv'), ...
        'Delimiter', ';');
else
    T_anova_all = table();
end

if ~isempty(all_posthoc_cond)
    T_posthoc_cond_all = vertcat(all_posthoc_cond{:});
    writetable(T_posthoc_cond_all, fullfile(globalOutDir, 'all_rois_posthoc_main_condition.csv'), ...
        'Delimiter', ';');
end

if ~isempty(all_posthoc_int)
    T_posthoc_int_all = vertcat(all_posthoc_int{:});
    writetable(T_posthoc_int_all, fullfile(globalOutDir, 'all_rois_posthoc_interaction.csv'), ...
        'Delimiter', ';');
end

% --- baseline t-tests ---
if ~isempty(all_baseline_rows)
    T_baseline_all = cell2table(all_baseline_rows, ...
        'VariableNames', {'ROI','Condition','Category','n','mean','sd','t','df','p_uncorr','p_bonf'});
    writetable(T_baseline_all, fullfile(globalOutDir, 'all_rois_baseline_tests.csv'), ...
        'Delimiter', ';');
else
    T_baseline_all = table();
end

% --- pairwise highest vs second-highest ---
if ~isempty(all_pair_rows)
    T_pair_all = cell2table(all_pair_rows, ...
        'VariableNames', {'ROI','Condition','Highest','SecondHighest','n','t','df','p_uncorr','p_bonf'});
    writetable(T_pair_all, fullfile(globalOutDir, 'all_rois_pairwise_highest_vs_second.csv'), ...
        'Delimiter', ';');
else
    T_pair_all = table();
end

%% ============================================================
%% PRINT ALL COMBINED TABLES TO COMMAND WINDOW
%% ============================================================
fprintf('\n============================================================\n');
fprintf('TABLE: Descriptives (M, SD per ROI x Condition x Category)\n');
fprintf('============================================================\n');
disp(T_desc_all);

fprintf('\n============================================================\n');
fprintf('TABLE: ANOVA Summary (all ROIs, Bonferroni N=%d per effect)\n', nAnovaCorr);
fprintf('============================================================\n');
disp(T_anova_all);

fprintf('\n============================================================\n');
fprintf('TABLE: Baseline t-tests vs Zero (Bonferroni N=%d, 3 ROIs x 5 categories)\n', nBaseCorr);
fprintf('============================================================\n');
disp(T_baseline_all);

fprintf('\n============================================================\n');
fprintf('TABLE: Pairwise Highest vs Second-Highest d'' (Bonferroni N=%d)\n', nPairCorr);
fprintf('============================================================\n');
disp(T_pair_all);

%% ============================================================
%% LATEX EXPORT: DESCRIPTIVES TABLE
%% Matches master_dprime_analysis_complete.m format:
%%   15 columns — ROI, Cat, n, Nat M/SD, Gat M/SD, Tex M/SD,
%%                N-G p, N-G p_bonf, N-T p, N-T p_bonf, G-T p, G-T p_bonf
%% ============================================================
if ~isempty(T_desc_all)

    latexFileDesc = fullfile(globalOutDir, sprintf('%s_all_rois_descriptives_table.tex', model_name));
    fid = fopen(latexFileDesc, 'w');

    fprintf(fid, '%% Auto-generated from MATLAB\n');
    fprintf(fid, '%% Condition post-hoc Bonferroni N=%d (3 pairs per ROI x category cell)\n\n', nDescPosthocCorr);
    fprintf(fid, '\\begin{table}[!ht]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\resizebox{\\textwidth}{!}{%%\n');
    fprintf(fid, '\\begin{tabular}{llc cc cc cc cccccc}\n');
    fprintf(fid, '\\toprule\n\n');

    % Row 1: top-level spanning headers
    fprintf(fid, '\\textbf{ROI} & \\textbf{Category} & \\textbf{n} ');
    fprintf(fid, '& \\multicolumn{6}{c}{\\textbf{Image Condition}} ');
    fprintf(fid, '& \\multicolumn{6}{c}{\\textbf{Post hoc}} \\\\\n\n');

    fprintf(fid, '\\cmidrule(lr){4-9} \\cmidrule(lr){10-15}\n\n');

    % Row 2: condition sub-headers and pair sub-headers
    fprintf(fid, ' &  &  ');
    fprintf(fid, '& \\multicolumn{2}{c}{\\textbf{Natural}} ');
    fprintf(fid, '& \\multicolumn{2}{c}{\\textbf{Gatys}} ');
    fprintf(fid, '& \\multicolumn{2}{c}{\\textbf{Texform}} ');
    fprintf(fid, '& \\multicolumn{2}{c}{\\textbf{N--G}} ');
    fprintf(fid, '& \\multicolumn{2}{c}{\\textbf{N--T}} ');
    fprintf(fid, '& \\multicolumn{2}{c}{\\textbf{G--T}} \\\\\n\n');

    fprintf(fid, ['\\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9} ' ...
                  '\\cmidrule(lr){10-11} \\cmidrule(lr){12-13} \\cmidrule(lr){14-15}\n\n']);

    % Row 3: M/SD and p/p_bonf sub-headers
    fprintf(fid, ' &  &  ');
    fprintf(fid, '& \\textbf{M} & \\textbf{SD} ');
    fprintf(fid, '& \\textbf{M} & \\textbf{SD} ');
    fprintf(fid, '& \\textbf{M} & \\textbf{SD} ');
    fprintf(fid, ['& \\textbf{p} & \\textbf{p\\textsubscript{Bonf=%d}} ' ...
                  '& \\textbf{p} & \\textbf{p\\textsubscript{Bonf=%d}} ' ...
                  '& \\textbf{p} & \\textbf{p\\textsubscript{Bonf=%d}} \\\\\n\n'], ...
                  nDescPosthocCorr, nDescPosthocCorr, nDescPosthocCorr);

    fprintf(fid, '\\midrule\n\n');

    roi_order = known_rois;
    cat_order = known_categories;

    for r = 1:numel(roi_order)
        roi_name = roi_order{r};
        roi_mask_tex = strcmp(T_desc_all.ROI, roi_name);
        T_roi = T_desc_all(roi_mask_tex, :);

        [~, idx] = ismember(T_roi.Category, cat_order);
        [~, sortIdx] = sort(idx);
        T_roi = T_roi(sortIdx, :);

        pretty_roi = localPrettyROI(roi_name);

        for i = 1:height(T_roi)
            if i == 1
                fprintf(fid, '\\multirow{5}{*}{\\textbf{%s}} ', pretty_roi);
            else
                fprintf(fid, ' ');
            end

            fprintf(fid, '& %s & %d & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %s & %s & %s & %s & %s \\\\\n', ...
                char(T_roi.Category{i}), ...
                T_roi.n_stimuli(i), ...
                T_roi.natural_mean(i), T_roi.natural_sd(i), ...
                T_roi.gatys_mean(i),   T_roi.gatys_sd(i), ...
                T_roi.texform_mean(i), T_roi.texform_sd(i), ...
                format_p_label(T_roi.p_NG(i)),      format_p_label(T_roi.p_NG_bonf(i)), ...
                format_p_label(T_roi.p_NT(i)),      format_p_label(T_roi.p_NT_bonf(i)), ...
                format_p_label(T_roi.p_GT(i)),      format_p_label(T_roi.p_GT_bonf(i)));
        end

        if r < numel(roi_order)
            fprintf(fid, '\n\\midrule\n\n');
        end
    end

    fprintf(fid, '\n\\bottomrule\n');
    fprintf(fid, '\\end{tabular}}\n');
    fprintf(fid, '\\end{table}\n');

    fclose(fid);
    fprintf('Saved LaTeX descriptives table: %s\n', latexFileDesc);
end

%% ============================================================
%% LATEX EXPORT: ANOVA TABLE
%% Matches master_dprime_analysis_complete.m format:
%%   Columns: ROI, Effect, F, df1, df2, p, p_bonf(=3), eta2
%% ============================================================
if ~isempty(T_anova_all)

    latexFileAnova = fullfile(globalOutDir, sprintf('%s_all_rois_anova_table.tex', model_name));
    fid = fopen(latexFileAnova, 'w');

    fprintf(fid, '%% Auto-generated from MATLAB\n');
    fprintf(fid, '%% Bonferroni correction across %d ROIs per effect type\n\n', nAnovaCorr);
    fprintf(fid, '\\begin{table}[!ht]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\begin{tabular}{llcccccc}\n');
    fprintf(fid, '\\toprule\n\n');

    fprintf(fid, ['\\textbf{ROI} & \\textbf{Effect} & \\textbf{F} & \\textbf{df1} ' ...
                  '& \\textbf{df2} & \\textbf{p} ' ...
                  '& \\textbf{p\\textsubscript{Bonf=%d}} & \\textbf{$\\eta^2$} \\\\\n\n'], ...
                  nAnovaCorr);
    fprintf(fid, '\\midrule\n\n');

    roi_order    = known_rois;
    effect_order = {'condition','category','condition:category'};

    for r = 1:numel(roi_order)
        roi_name = roi_order{r};
        roi_mask_tex = strcmp(T_anova_all.ROI, roi_name);
        T_roi = T_anova_all(roi_mask_tex, :);

        [~, idx] = ismember(T_roi.Effect, effect_order);
        [~, sortIdx] = sort(idx);
        T_roi = T_roi(sortIdx, :);

        pretty_roi = localPrettyROI(roi_name);

        for i = 1:height(T_roi)
            if i == 1
                fprintf(fid, '\\multirow{3}{*}{\\textbf{%s}} ', pretty_roi);
            else
                fprintf(fid, ' ');
            end

            fprintf(fid, '& %s & %.2f & %d & %d & %s & %s & %.3f \\\\\n', ...
                char(T_roi.Effect{i}), ...
                T_roi.F(i), ...
                T_roi.df1(i), ...
                T_roi.df2(i), ...
                format_p_label_plain(T_roi.p(i)), ...
                format_p_label_plain(T_roi.p_bonf(i)), ...
                T_roi.eta2(i));
        end

        if r < numel(roi_order)
            fprintf(fid, '\n\\midrule\n\n');
        end
    end

    fprintf(fid, '\n\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');

    fclose(fid);
    fprintf('Saved LaTeX ANOVA table: %s\n', latexFileAnova);
end

%% ============================================================
%% LATEX EXPORT: BASELINE T-TEST TABLE
%% Matches master_dprime_analysis_complete.m: fmri_all_rois_baseline_table.tex
%% Columns: ROI, Condition, Category, n, M, SD, t, df, p, p_bonf(=15)
%% ============================================================
if ~isempty(T_baseline_all)

    latexFileBase = fullfile(globalOutDir, sprintf('%s_all_rois_baseline_table.tex', model_name));
    fid = fopen(latexFileBase, 'w');

    fprintf(fid, '%% Auto-generated from MATLAB\n');
    fprintf(fid, '%% One-sample t-tests against zero per ROI x Condition x Category\n');
    fprintf(fid, '%% Bonferroni correction across %d comparisons (3 ROIs x 5 categories) per condition\n\n', nBaseCorr);
    fprintf(fid, '\\begin{table}[!ht]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\resizebox{\\textwidth}{!}{%%\n');
    fprintf(fid, '\\begin{tabular}{lllccccccc}\n');
    fprintf(fid, '\\toprule\n\n');

    fprintf(fid, ['\\textbf{ROI} & \\textbf{Condition} & \\textbf{Category} ' ...
                  '& \\textbf{n} & \\textbf{M} & \\textbf{SD} ' ...
                  '& \\textbf{t} & \\textbf{df} ' ...
                  '& \\textbf{p} ' ...
                  '& \\textbf{p\\textsubscript{Bonf=%d}} \\\\\n\n'], nBaseCorr);
    fprintf(fid, '\\midrule\n\n');

    condLabelMap = containers.Map( ...
        {'natural','gatys','texform'}, ...
        {'Natural','Gatys','Texform'});

    for r = 1:numel(known_rois)
        roi_name  = known_rois{r};
        nRowsROI  = n_cond * numel(known_categories);

        first_roi_row = true;
        for ci = 1:n_cond
            cond_name  = known_conditions{ci};
            nRowsCond  = numel(known_categories);

            first_cond_row = true;
            for j = 1:numel(known_categories)
                this_cat = known_categories{j};
                rowMask  = strcmp(T_baseline_all.ROI, roi_name) ...
                         & strcmp(T_baseline_all.Condition, cond_name) ...
                         & strcmp(T_baseline_all.Category, this_cat);
                row = T_baseline_all(rowMask, :);

                if isempty(row); continue; end

                if first_roi_row
                    fprintf(fid, '\\multirow{%d}{*}{\\textbf{%s}} ', nRowsROI, localPrettyROI(roi_name));
                    first_roi_row = false;
                else
                    fprintf(fid, ' ');
                end

                if first_cond_row
                    fprintf(fid, '& \\multirow{%d}{*}{%s} ', nRowsCond, condLabelMap(cond_name));
                    first_cond_row = false;
                else
                    fprintf(fid, '&  ');
                end

                fprintf(fid, '& %s & %d & %.2f & %.2f & %.2f & %d & %s & %s \\\\\n', ...
                    this_cat, row.n, row.mean, row.sd, row.t, row.df, ...
                    format_p_label_plain(row.p_uncorr), ...
                    format_p_label_plain(row.p_bonf));
            end

            if ci < n_cond
                fprintf(fid, '\\cmidrule(l){2-10}\n');
            end
        end

        if r < numel(known_rois)
            fprintf(fid, '\n\\midrule\n\n');
        end
    end

    fprintf(fid, '\n\\bottomrule\n');
    fprintf(fid, '\\end{tabular}}\n');
    fprintf(fid, ['\\caption{One-sample $t$-tests of the selectivity index (d$''$) against zero ' ...
                  'for each ROI, image condition, and stimulus category. ' ...
                  'p-values are reported uncorrected and Bonferroni-corrected across ' ...
                  '%d comparisons (3 ROIs $\\times$ 5 categories) per condition.}\n'], nBaseCorr);
    fprintf(fid, '\\end{table}\n');

    fclose(fid);
    fprintf('Saved LaTeX baseline table: %s\n', latexFileBase);
end

%% ============================================================
%% LATEX EXPORT: PAIRWISE HIGHEST VS SECOND-HIGHEST TABLE
%% ============================================================
if ~isempty(T_pair_all)

    latexFilePair = fullfile(globalOutDir, sprintf('%s_all_rois_pairwise_table.tex', model_name));
    fid = fopen(latexFilePair, 'w');

    fprintf(fid, '%% Auto-generated from MATLAB\n');
    fprintf(fid, '%% Pairwise test: highest d'' category vs second-highest, per ROI x Condition\n');
    fprintf(fid, '%% Independent-samples t-test (stimuli are independent across categories)\n');
    fprintf(fid, '%% Bonferroni correction across %d comparisons (3 conditions per ROI)\n\n', nPairCorr);
    fprintf(fid, '\\begin{table}[!ht]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\begin{tabular}{llllccccc}\n');
    fprintf(fid, '\\toprule\n\n');

    fprintf(fid, ['\\textbf{ROI} & \\textbf{Condition} & \\textbf{Highest} ' ...
                  '& \\textbf{2nd Highest} & \\textbf{n} & \\textbf{t} ' ...
                  '& \\textbf{df} & \\textbf{p} ' ...
                  '& \\textbf{p\\textsubscript{Bonf=%d}} \\\\\n\n'], nPairCorr);
    fprintf(fid, '\\midrule\n\n');

    condLabelMap = containers.Map( ...
        {'natural','gatys','texform'}, ...
        {'Natural','Gatys','Texform'});

    for r = 1:numel(known_rois)
        roi_name = known_rois{r};

        for ci = 1:n_cond
            cond_name = known_conditions{ci};
            rowMask   = strcmp(T_pair_all.ROI, roi_name) & strcmp(T_pair_all.Condition, cond_name);
            row       = T_pair_all(rowMask, :);

            if isempty(row); continue; end

            if ci == 1
                fprintf(fid, '\\multirow{%d}{*}{\\textbf{%s}} ', n_cond, localPrettyROI(roi_name));
            else
                fprintf(fid, ' ');
            end

            fprintf(fid, '& %s & %s & %s & %d & %.2f & %d & %s & %s \\\\\n', ...
                condLabelMap(cond_name), ...
                char(row.Highest), ...
                char(row.SecondHighest), ...
                row.n, row.t, row.df, ...
                format_p_label_plain(row.p_uncorr), ...
                format_p_label_plain(row.p_bonf));
        end

        if r < numel(known_rois)
            fprintf(fid, '\n\\midrule\n\n');
        end
    end

    fprintf(fid, '\n\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
   
    fprintf(fid, '\\end{table}\n');

    fclose(fid);
    fprintf('Saved LaTeX pairwise table: %s\n', latexFilePair);
end

fprintf('\n============================================================\n');
fprintf('Done.\n');
fprintf('Combined tables saved in:\n%s\n', globalOutDir);
fprintf('Model label used for LaTeX filenames: %s\n', model_name);
fprintf('============================================================\n');

%% ============================================================
%% LOCAL HELPER FUNCTIONS
%% ============================================================
function pretty_roi = localPrettyROI(roi_name)
    switch roi_name
        case 'LOTC_body_left'
            pretty_roi = 'LOTC Body Left';
        case 'LOTC_hand_left'
            pretty_roi = 'LOTC Hand Left';
        case 'LOTC_tool_left'
            pretty_roi = 'LOTC Tool Left';
        otherwise
            pretty_roi = strrep(roi_name, '_', ' ');
    end
end

function p_str = format_p_label(p)

    if isnan(p)
        p_str = 'n/a';
    elseif p < 0.001
        p_str = '$p {<} .001$';
    else
        p_rounded = round(p, 3);
        p_str = sprintf('$p = .%03d$', round(p_rounded * 1000));
    end
end

function p_str = format_p_label_plain(p)
  
    if isnan(p)
        p_str = 'n/a';
    elseif p < 0.001
        p_str = '{$<$}\,.001';
    else
        p_rounded = round(p, 3);
        p_str = sprintf('.%03d', round(p_rounded * 1000));
    end
end
