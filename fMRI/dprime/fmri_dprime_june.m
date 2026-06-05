
clear; clc; close all;

%% ============================================================
%% 1. PATHS & SETTINGS
%% ============================================================
baseDir   = '/Users/lunameidoering/Desktop/CODE/fMRI/multivariate/';
outputDir = '/Users/lunameidoering/Desktop/CODE/fMRI/dprime_final_results_0306/';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

subjects   = {'sub01', 'sub02', 'sub03', 'sub04', 'sub05','sub06', 'sub07', 'sub08', 'sub09', 'sub10', 'sub11', 'sub12', 'sub14'};
conditions = {'original','gatys','texform'};
rois       = {'LOTC_body','LOTC_hand','LOTC_tool'};
categories = {'body','hand','tool','manip','nonmanip'};
catNames   = {'Bodies', 'Hands', 'Tools', 'Manip', 'Nman'};

catColors  = [1 0.5 0; 0.3 0.6 0.3; 0 0 0.5; 0.3 0.75 0.93; 1 0 1];
yLims      = [-2, 4]; % Fixed standard Y-axis limits for all plots

% Bonferroni corection values
nBaseCorr         = 15;  % one-sample vs-zero: 3 ROIs x 5 categories (per condition)
nPairCorr         = 3;   % pairwise highest vs second-highest: 3 conditions (per ROI)
nAnovaCorr        = 3;   
nDescPosthocCorr  = 3;   

condFolders = containers.Map({'original','gatys','texform'}, {'ds_original','ds_gatys','ds_texform'});

% indexing logic
% Rows in S.ds.samples are trials, columns are voxels
idx_gatys_texform = {[1,6,11,16], [2,7,12,17], [3,8,13,18], [4,9,14,19], [5,10,15,20]};
idx_orig_default  = {[1,6,11], [2,7,12], [3,8,13], [4,9,14], [5,10,15]};
idx_orig_short    = {[1,6], [2,7], [3,8], [4,9], [5,10]};
subjectsWithFewerorig = {'sub02','sub04','sub12'};

%% ============================================================
%% 2. DATA EXTRACTION & D-PRIME CALCULATION
%% ============================================================
nCond = numel(conditions);
nROI  = numel(rois);
nCat  = numel(categories);
nSubj = numel(subjects);
nMeas = nCond * nROI * nCat;

dataMat_dprime   = nan(nSubj, nMeas);
measureNames     = cell(1, nMeas);
withinCondFactor = cell(nMeas,1);
withinROIFactor  = cell(nMeas,1);
withinCatFactor  = cell(nMeas,1);

% for csvs
subjectCsvDir = fullfile(outputDir, 'subject_csvs_for_python');
if ~exist(subjectCsvDir, 'dir'), mkdir(subjectCsvDir); end 


%% COLLECTORS FOR COMBINED TABLES
all_desc_rows = {};
all_anova_rows = {};
all_desc_posthoc_rows = {};
all_pair_rows      = {};   % pairwise: target vs highest (or highest vs second if target=highest)
all_baseline_rows  = {};   % per-category t-tests against zero

for r = 1:nROI
    skip_rows = {}; % log skipped data
    thisROI = rois{r};

    % Target category index: the category whose selectivity this ROI is named for
    switch thisROI
        case 'LOTC_body', targetCatIdx = 1;   % 'body'
        case 'LOTC_hand', targetCatIdx = 2;   % 'hand'
        case 'LOTC_tool', targetCatIdx = 3;   % 'tool'
        otherwise,        targetCatIdx = 1;
    end

    for c = 1:nCond
        thisCond = conditions{c};
        roi_cond_dprime_vals = cell(1, nCat);

        for s = 1:nSubj
            subj = subjects{s};
            condPath = fullfile(baseDir, condFolders(thisCond));
            pattern  = sprintf('%s*%s*_ds.mat', subj, thisROI);
            f = dir(fullfile(condPath, pattern));

            if isempty(f)
                msg = sprintf('No file found for subject %s, ROI %s, condition %s', subj, thisROI, thisCond);
                warning(msg);
            
                skip_rows(end+1,:) = {subj, thisROI, thisCond, 'missing_file', NaN, NaN, msg}; 
                continue;

        
            elseif numel(f) > 1
                warning('Multiple files found for subject %s, ROI %s, condition %s. Using first match: %s', ...
                    subj, thisROI, thisCond, f(1).name);
            end

            S = load(fullfile(condPath, f(1).name));
            if ~isfield(S, 'ds') || ~isfield(S.ds, 'samples')
                msg = sprintf('Missing S.ds.samples in file: %s', fullfile(condPath, f(1).name));
                warning(msg);
            
                skip_rows(end+1,:) = {subj, thisROI, thisCond, 'missing_ds_samples', NaN, NaN, msg}; 
                continue;
            end

            trial_by_voxel = S.ds.samples;   

            if strcmp(thisCond, 'original')
                if ismember(subj, subjectsWithFewerorig)
                    currIdx = idx_orig_short;
                else
                    currIdx = idx_orig_default;
                end
            else
                currIdx = idx_gatys_texform;
            end

            % Basic consistency checks
            nTrials = size(trial_by_voxel, 1);
            maxIdxNeeded = max([currIdx{:}]);
            if nTrials < maxIdxNeeded
                msg = sprintf(['Not enough trials in file for subject %s, ROI %s, condition %s. ' ...
                               'Expected at least %d trials, found %d.'], ...
                               subj, thisROI, thisCond, maxIdxNeeded, nTrials);
            
                warning(msg);
            
                skip_rows(end+1,:) = {subj, thisROI, thisCond, 'not_enough_trials', maxIdxNeeded, nTrials, msg}; %#ok<SAGROW>
                continue;
            end

            % ROI mean activation per trial 
            all_means = mean(trial_by_voxel, 2, 'omitnan');   % nTrials x 1

          
            % d' = (mu_cat - mu_noncat) / sqrt(0.5 * (var_cat + var_noncat)
            for g = 1:nCat
                data_cat = all_means(currIdx{g});

                otherCats = setdiff(1:nCat, g);
                otherIdx  = [currIdx{otherCats}];
                data_noncat = all_means(otherIdx);

                mu_cat    = mean(data_cat, 'omitnan');
                mu_noncat = mean(data_noncat, 'omitnan');

                var_cat    = var(data_cat, 'omitnan');
                var_noncat = var(data_noncat, 'omitnan');

                sigma_pooled = sqrt(0.5 * (var_cat + var_noncat));

                if isnan(sigma_pooled) || sigma_pooled == 0
                    d_val = NaN;
                else
                    d_val = (mu_cat - mu_noncat) / sigma_pooled;
                end

                roi_cond_dprime_vals{g}(s,1) = d_val;

                colIdx = (c-1)*nROI*nCat + (r-1)*nCat + g;
                dataMat_dprime(s, colIdx) = d_val;

                if s == 1
                    measureNames{colIdx}     = matlab.lang.makeValidName(sprintf('%s_%s_%s', thisCond, thisROI, categories{g}));
                    withinCondFactor{colIdx} = thisCond;
                    withinROIFactor{colIdx}  = thisROI;
                    withinCatFactor{colIdx}  = categories{g};
                end
            end
        end


        
        % save per subject d prime averages
        for s = 1:nSubj
            subj = subjects{s};

            subj_vals = nan(1, nCat);
    
            for g = 1:nCat
                if numel(roi_cond_dprime_vals{g}) >= s
                    subj_vals(g) = roi_cond_dprime_vals{g}(s);
                else
                    subj_vals(g) = NaN;
                end
            end

            % Convert ROI name to short label 
            switch thisROI
                case 'LOTC_body'
                    roi_short = 'body';
                case 'LOTC_hand'
                    roi_short = 'hand';
                case 'LOTC_tool'
                    roi_short = 'tool';
                otherwise
                    roi_short = strrep(thisROI, 'LOTC_', '');
            end

          
            T_subj = array2table(subj_vals, ...
                'VariableNames', {'bodies','hands','tools','mani','nman'});

            % Save
            outName = sprintf('%s_%s_%s_activations.csv', subj, thisCond, roi_short);
            writetable(T_subj, fullfile(subjectCsvDir, outName));
        end

        
        % Stats and barplots for current ROI/Condition
        m_mean = cellfun(@(x) mean(x, 'omitnan'), roi_cond_dprime_vals);
        m_sem  = cellfun(@(x) std(x, 'omitnan') / sqrt(sum(~isnan(x))), roi_cond_dprime_vals);
        m_sd   = cellfun(@(x) std(x, 'omitnan'), roi_cond_dprime_vals);
        n_vals = cellfun(@(x) sum(~isnan(x)), roi_cond_dprime_vals);

        test_vs0  = test_against_zero(roi_cond_dprime_vals, nCat, nBaseCorr);
        pair_test = test_highest_vs_second(roi_cond_dprime_vals, m_mean, catNames, targetCatIdx);
        pair_test.p_bonf = min(pair_test.p * nPairCorr, 1);   % corrected for nPairCorr tests

        file_tag = sprintf('%s_%s', thisROI, thisCond);
        save_metric_csvs(outputDir, file_tag, catNames, n_vals, m_mean, m_sem, m_sd, test_vs0, pair_test);

        % collect descriptives for combined table
        for g = 1:nCat
            all_desc_rows(end+1, :) = { ...
                thisROI, ...
                thisCond, ...
                categories{g}, ...
                n_vals(g), ...
                m_mean(g), ...
                m_sem(g), ...
                m_sd(g), ...
                test_vs0.p_bonf(g) ...
            }; 
        end

        for g = 1:nCat
            all_baseline_rows(end+1, :) = { ...
                thisROI, thisCond, categories{g}, ...
                n_vals(g), m_mean(g), m_sd(g), ...
                test_vs0.t(g), test_vs0.df(g), ...
                test_vs0.p(g), test_vs0.p_bonf(g) ...
            }; 
        end

        % collect pairwise rows
        all_pair_rows(end+1, :) = { ...
            thisROI, thisCond, ...
            pair_test.cat1, pair_test.cat2, pair_test.isTargetHighest, ...
            pair_test.n, pair_test.t, pair_test.df, pair_test.p, pair_test.p_bonf ...
        }; 

        outFile = fullfile(outputDir, [file_tag '_barplot.png']);
        plot_metric_barplot(m_mean, m_sem, test_vs0.p_bonf, pair_test, catNames, catColors, ...
            strrep(file_tag,'_',' '), "Selectivity Index (d')", outFile, yLims);
    end
end

%% ============================================================
%% SAVE SKIP LOG
%% ============================================================
if ~isempty(skip_rows)
    skipTable = cell2table(skip_rows, ...
        'VariableNames', {'Subject','ROI','Condition','Reason','ExpectedTrials','FoundTrials','Message'});

    writetable(skipTable, fullfile(outputDir, 'skipped_data_log.csv'));

    fprintf('\n============================================================\n');
    fprintf('Skipped data summary:\n');
    fprintf('============================================================\n');
    disp(skipTable);
else
    fprintf('\nNo data were skipped.\n');
end

%% ============================================================
%% 3. ANOVA, POST-HOCS, & INTERACTION PLOTS
%% ============================================================



wideTable = array2table(dataMat_dprime, 'VariableNames', measureNames);
wideTable = addvars(wideTable, subjects', 'Before', 1, 'NewVariableNames', 'Subject');

withinDesign = table( ...
    categorical(withinCondFactor, conditions, 'Ordinal', true), ...
    categorical(withinROIFactor, rois, 'Ordinal', true), ...
    categorical(withinCatFactor, categories, 'Ordinal', true), ...
    'VariableNames', {'Condition','ROI','Category'});

for r = 1:nROI
    thisROI = rois{r};

    fprintf('\n=======================================================\n');
    fprintf('=== RM-ANOVA RESULTS FOR ROI: %s ===\n', thisROI);
    fprintf('=======================================================\n');

    roiMask = strcmp(cellstr(withinDesign.ROI), thisROI);
    thisMeasureNames = measureNames(roiMask);
    thisWithin = withinDesign(roiMask, :);
    thisWide = wideTable(:, ['Subject', thisMeasureNames]);


    %% ------------------------------------------------------------
    %% REMOVE SUBJECTS WITH MISSING DATA FOR THIS ROI-WISE ANOVA
    %% ------------------------------------------------------------
    dataOnlyForROI = thisWide{:, thisMeasureNames};
    completeMask = all(~isnan(dataOnlyForROI), 2);
    
    if any(~completeMask)
        skippedSubjects = thisWide.Subject(~completeMask);
    
        warning('Skipping %d subject(s) from ANOVA for ROI %s due to incomplete data.', ...
            sum(~completeMask), thisROI);
    
        disp(table(skippedSubjects, ...
            'VariableNames', {'SubjectsSkippedFromANOVA'}));
    end
    
    thisWide = thisWide(completeMask, :);
    
    if height(thisWide) < 2
        warning('Not enough complete subjects for ANOVA in ROI %s. Skipping ANOVA.', thisROI);
        continue;
    end

    measRange = sprintf('%s-%s ~ 1', thisMeasureNames{1}, thisMeasureNames{end});
    rm_roi = fitrm(thisWide, measRange, 'WithinDesign', thisWithin);
    ranova_roi = ranova(rm_roi, 'WithinModel', 'Condition*Category');

    % Print ANOVA to command window and save
    disp(ranova_roi);
    writetable(ranova_roi, fullfile(outputDir, sprintf('ranova_dprime_%s.csv', thisROI)));


    % collect ANOVA rows for combined export
    ranova_names = ranova_roi.Properties.RowNames;

    idx_cond     = strcmp(ranova_names, '(Intercept):Condition');
    idx_err_cond = strcmp(ranova_names, 'Error(Condition)');

    idx_cat      = strcmp(ranova_names, '(Intercept):Category');
    idx_err_cat  = strcmp(ranova_names, 'Error(Category)');

    idx_int      = strcmp(ranova_names, '(Intercept):Condition:Category');
    idx_err_int  = strcmp(ranova_names, 'Error(Condition:Category)');

    if any(idx_cond) && any(idx_err_cond)
        all_anova_rows(end+1, :) = { ...
            thisROI, 'condition', ...
            ranova_roi.F(idx_cond), ...
            ranova_roi.DF(idx_cond), ...
            ranova_roi.DF(idx_err_cond), ...
            ranova_roi.pValueGG(idx_cond) ...
        }; %#ok<SAGROW>
    end

    if any(idx_cat) && any(idx_err_cat)
        all_anova_rows(end+1, :) = { ...
            thisROI, 'category', ...
            ranova_roi.F(idx_cat), ...
            ranova_roi.DF(idx_cat), ...
            ranova_roi.DF(idx_err_cat), ...
            ranova_roi.pValueGG(idx_cat) ...
        }; %#ok<SAGROW>
    end

    if any(idx_int) && any(idx_err_int)
        all_anova_rows(end+1, :) = { ...
            thisROI, 'condition:category', ...
            ranova_roi.F(idx_int), ...
            ranova_roi.DF(idx_int), ...
            ranova_roi.DF(idx_err_int), ...
            ranova_roi.pValueGG(idx_int) ...
        }; %#ok<SAGROW>
    end


    % Reshape for plots and post-hocs
    dataOnly = thisWide{:, thisMeasureNames};      % nSubj x (nCond*nCat)
    nSubjROI = size(dataOnly, 1);
    Yroi = reshape(dataOnly', [nCat, nCond, nSubjROI]);
    Yroi = permute(Yroi, [3 2 1]);                 % nSubjROI x nCond x nCat

    %% --- DESCRIPTIVE TABLE POST-HOCS (paired condition tests within each category) ---
    for g = 1:nCat
        v_orig = Yroi(:, 1, g);
        v_gat  = Yroi(:, 2, g);
        v_tex  = Yroi(:, 3, g);

        % remove NaN pairs jointly
        mask_og = ~isnan(v_orig) & ~isnan(v_gat);
        mask_ot = ~isnan(v_orig) & ~isnan(v_tex);
        mask_gt = ~isnan(v_gat)  & ~isnan(v_tex);

        if sum(mask_og) >= 2
            [~, p_og] = ttest(v_orig(mask_og), v_gat(mask_og));
        else
            p_og = NaN;
        end

        if sum(mask_ot) >= 2
            [~, p_ot] = ttest(v_orig(mask_ot), v_tex(mask_ot));
        else
            p_ot = NaN;
        end

        if sum(mask_gt) >= 2
            [~, p_gt] = ttest(v_gat(mask_gt), v_tex(mask_gt));
        else
            p_gt = NaN;
        end

        all_desc_posthoc_rows(end+1, :) = { ...
            thisROI, categories{g}, ...
            p_og, min(p_og * nDescPosthocCorr, 1), ...
            p_ot, min(p_ot * nDescPosthocCorr, 1), ...
            p_gt, min(p_gt * nDescPosthocCorr, 1) ...
        }; %#ok<SAGROW>
    end

    %% --- CALCULATE & PRINT POST-HOCS ---
    fprintf('\n--- Interaction Post-Hocs (Condition Pairs per Category) ---\n');
    intRows = {};
    pairs = nchoosek(1:nCond, 2);
    nPairs = size(pairs, 1);

    for g = 1:nCat
        for i = 1:nPairs
            a = pairs(i,1);
            b = pairs(i,2);

            [~, p, ~, stats] = ttest(Yroi(:,a,g), Yroi(:,b,g));

            if p < 0.001
                sig_label = '***';
            elseif p < 0.01
                sig_label = '**';
            elseif p < 0.05
                sig_label = '*';
            else
                sig_label = 'n.s.';
            end

            intRows(end+1,:) = {catNames{g}, conditions{a}, conditions{b}, stats.tstat, p, sig_label}; %#ok<AGROW>
        end
    end

    T_int = cell2table(intRows, 'VariableNames', {'Category','C1','C2','t','p_raw','sig'});
    T_int.p_bonf = min(T_int.p_raw * (nCat * nPairs), 1);
    disp(T_int);
    writetable(T_int, fullfile(outputDir, sprintf('table_%s_posthoc.csv', thisROI)));

    %% --- PLOTS ---
    intMean = squeeze(mean(Yroi, 1, 'omitnan'));
    intSEM  = squeeze(std(Yroi, 0, 1, 'omitnan')) / sqrt(nSubjROI);

    catSubjData  = squeeze(mean(Yroi, 2, 'omitnan')); % across conditions
    condSubjData = squeeze(mean(Yroi, 3, 'omitnan')); % across categories

    % --- 1. Interaction Plot ---
    fig1 = figure('Position', [100 100 700 500], 'Visible', 'off');
    hold on;

    for g = 1:nCat
        errorbar(1:nCond, intMean(:,g), intSEM(:,g), '-o', ...
            'Color', catColors(g,:), ...
            'LineWidth', 2, ...
            'MarkerSize', 7, ...
            'MarkerFaceColor', catColors(g,:), ...
            'CapSize', 8);
    end

    set(gca, 'XTick', 1:nCond, 'XTickLabel', conditions, 'FontSize', 20);


    box off;

    exportgraphics(fig1, fullfile(outputDir, sprintf('plot_%s_interaction.png', thisROI)), 'Resolution', 300);

    % --- 2. Condition Main Effect Plot ---
    fig2 = figure('Position', [100 100 600 450], 'Visible', 'off');
    hold on;

    plot(1:nCond, condSubjData', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);

    errorbar(1:nCond, mean(condSubjData, 1, 'omitnan'), std(condSubjData, 0, 1, 'omitnan')/sqrt(nSubjROI), ...
        '-ko', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'CapSize', 10);

    set(gca, 'XTick', 1:nCond, 'XTickLabel', conditions, 'FontSize', 34);
    %xlabel('Image Condition', 'FontSize', 12, 'FontWeight', 'bold');
    %ylabel("Average d'", 'FontSize', 12, 'FontWeight', 'bold');
    %title(sprintf('%s: Condition Main Effect', strrep(thisROI, '_', ' ')), 'Interpreter', 'none');

    grid on;
    box off;

    exportgraphics(fig2, fullfile(outputDir, sprintf('plot_%s_main_cond.png', thisROI)), 'Resolution', 300);

    close(fig1);
    close(fig2);
end

save(fullfile(outputDir, 'rm_anova_tables_dprime.mat'), 'wideTable', 'withinDesign', 'measureNames');
fprintf('\nDone. Results in: %s\n', outputDir);


%% ============================================================
%% BUILD COMBINED TABLES
%% ============================================================
T_desc_all = cell2table(all_desc_rows, ...
    'VariableNames', {'ROI','Condition','Category','n','mean','sem','sd','p_bonf'});

writetable(T_desc_all, fullfile(outputDir, 'all_rois_descriptives.csv'));

T_desc_posthoc_all = cell2table(all_desc_posthoc_rows, ...
    'VariableNames', {'ROI','Category', ...
        'p_OG','p_OG_bonf', ...
        'p_OT','p_OT_bonf', ...
        'p_GT','p_GT_bonf'});

writetable(T_desc_posthoc_all, fullfile(outputDir, 'all_rois_descriptive_posthocs.csv'));

T_anova_all = cell2table(all_anova_rows, ...
    'VariableNames', {'ROI','Effect','F','df1','df2','p'});

T_anova_all.p_bonf = min(T_anova_all.p * nAnovaCorr, 1);

writetable(T_anova_all, fullfile(outputDir, 'all_rois_anova.csv'));

% --- COMBINED PAIRWISE TABLE ---

T_pair_all = cell2table(all_pair_rows, ...
    'VariableNames', {'ROI','Condition','Cat1','Cat2','TargetIsHighest','n','t','df','p_uncorr','p_bonf'});
writetable(T_pair_all, fullfile(outputDir, 'all_rois_pairwise_target_vs_highest.csv'));

T_baseline_all = cell2table(all_baseline_rows, ...
    'VariableNames', {'ROI','Condition','Category','n','mean','sd','t','df','p_uncorr','p_bonf'});
writetable(T_baseline_all, fullfile(outputDir, 'all_rois_baseline_tests.csv'));

%% --- PRINT ALL COMBINED TABLES TO COMMAND WINDOW ---
fprintf('\n============================================================\n');
fprintf('TABLE: Descriptives (mean, SD per ROI x Condition x Category)\n');
fprintf('============================================================\n');
disp(T_desc_all);

fprintf('\n============================================================\n');
fprintf('TABLE: Descriptive Post-Hocs (condition pair p-values per ROI x Category)\n');
fprintf('============================================================\n');
disp(T_desc_posthoc_all);

fprintf('\n============================================================\n');
fprintf('TABLE: RM-ANOVA Summary (all ROIs)\n');
fprintf('============================================================\n');
disp(T_anova_all);

fprintf('\n============================================================\n');
fprintf('TABLE: Pairwise target vs highest d'' (Bonferroni N=%d)\n', nPairCorr);
fprintf('  Cat1 = highest-d'' category; Cat2 = target ROI category\n');
fprintf('  (when target IS highest: Cat2 = second-highest instead)\n');
fprintf('============================================================\n');
disp(T_pair_all);

fprintf('\n============================================================\n');
fprintf('TABLE: Baseline t-tests vs Zero (Bonferroni N=%d, across 3 ROIs x 5 categories)\n', nBaseCorr);
fprintf('============================================================\n');
disp(T_baseline_all);

%% ============================================================
%% LATEX EXPORT: DESCRIPTIVES TABLE
%% ============================================================
if ~isempty(T_desc_all)

    latexFileDesc = fullfile(outputDir, 'fmri_all_rois_descriptives_table.tex');
    fid = fopen(latexFileDesc, 'w');

    fprintf(fid, '%% Auto-generated from MATLAB\n');
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

    roi_order = rois;
    cat_order = categories;

    for r = 1:numel(roi_order)
        roi_name = roi_order{r};
        T_roi = T_desc_all(strcmp(T_desc_all.ROI, roi_name), :);
        T_post = T_desc_posthoc_all(strcmp(T_desc_posthoc_all.ROI, roi_name), :);

        for i = 1:numel(cat_order)
            this_cat = cat_order{i};

            row_orig = T_roi(strcmp(T_roi.Category, this_cat) & strcmp(T_roi.Condition, 'original'), :);
            row_gat  = T_roi(strcmp(T_roi.Category, this_cat) & strcmp(T_roi.Condition, 'gatys'), :);
            row_tex  = T_roi(strcmp(T_roi.Category, this_cat) & strcmp(T_roi.Condition, 'texform'), :);
            row_post = T_post(strcmp(T_post.Category, this_cat), :);

            if i == 1
                fprintf(fid, '\\multirow{5}{*}{\\textbf{%s}} ', localPrettyROI(roi_name));
            else
                fprintf(fid, ' ');
            end

            fprintf(fid, '& %s & %d & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %s & %s & %s & %s & %s \\\\\n', ...
                localPrettyCat(this_cat), ...
                row_orig.n, ...
                row_orig.mean, row_orig.sd, ...
                row_gat.mean,  row_gat.sd, ...
                row_tex.mean,  row_tex.sd, ...
                format_p_label(row_post.p_OG),      format_p_label(row_post.p_OG_bonf), ...
                format_p_label(row_post.p_OT),      format_p_label(row_post.p_OT_bonf), ...
                format_p_label(row_post.p_GT),      format_p_label(row_post.p_GT_bonf));
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
%% ============================================================
if ~isempty(T_anova_all)

    latexFileAnova = fullfile(outputDir, 'fmri_all_rois_anova_table.tex');
    fid = fopen(latexFileAnova, 'w');

    fprintf(fid, '%% Auto-generated from MATLAB\n');
    fprintf(fid, '%% Bonferroni correction across %d ROIs per effect type\n\n', nAnovaCorr);
    fprintf(fid, '\\begin{table}[!ht]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\begin{tabular}{llccccc}\n');
    fprintf(fid, '\\toprule\n\n');

    fprintf(fid, ['\\textbf{ROI} & \\textbf{Effect} & \\textbf{F} & \\textbf{df1} ' ...
                  '& \\textbf{df2} & \\textbf{p} & \\textbf{p\\textsubscript{Bonf=%d}} \\\\\n\n'], ...
                  nAnovaCorr);
    fprintf(fid, '\\midrule\n\n');

    roi_order = rois;
    effect_order = {'condition','category','condition:category'};

    for r = 1:numel(roi_order)
        roi_name = roi_order{r};
        T_roi = T_anova_all(strcmp(T_anova_all.ROI, roi_name), :);

        [~, idx] = ismember(T_roi.Effect, effect_order);
        [~, sortIdx] = sort(idx);
        T_roi = T_roi(sortIdx, :);

        for i = 1:height(T_roi)
            if i == 1
                fprintf(fid, '\\multirow{3}{*}{\\textbf{%s}} ', localPrettyROI(roi_name));
            else
                fprintf(fid, ' ');
            end

            fprintf(fid, '& %s & %.2f & %d & %d & %s & %s \\\\\n', ...
                char(T_roi.Effect{i}), ...
                T_roi.F(i), ...
                T_roi.df1(i), ...
                T_roi.df2(i), ...
                format_p_label_plain(T_roi.p(i)), ...
                format_p_label_plain(T_roi.p_bonf(i)));
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
%% LATEX EXPORT: PAIRWISE TARGET VS HIGHEST TABLE
%% ============================================================
% Cat1 = highest-d' category; Cat2 = target ROI category
% When target IS the highest: Cat2 = second-highest (flagged by TargetIsHighest = true)
% Bonferroni N = nPairCorr (= 3, one per condition within an ROI)
if ~isempty(T_pair_all)

    latexFilePair = fullfile(outputDir, 'fmri_all_rois_pairwise_table.tex');
    fid = fopen(latexFilePair, 'w');

    fprintf(fid, '%% Auto-generated from MATLAB\n');
    fprintf(fid, '%% Pairwise test: highest d'' category vs target ROI category\n');
    fprintf(fid, '%% When target IS highest: highest vs second-highest (TargetIsHighest = true)\n');
    fprintf(fid, '%% Bonferroni correction across %d comparisons (3 conditions per ROI)\n\n', nPairCorr);
    fprintf(fid, '\\begin{table}[!ht]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\begin{tabular}{lllllccccc}\n');
    fprintf(fid, '\\toprule\n\n');
    fprintf(fid, ['\\textbf{ROI} & \\textbf{Condition} & \\textbf{Cat 1} ' ...
                  '& \\textbf{Cat 2} & \\textbf{Target = Highest} & \\textbf{n} & \\textbf{t} ' ...
                  '& \\textbf{df} & \\textbf{p} ' ...
                  '& \\textbf{p\\textsubscript{Bonf=%d}} \\\\\n\n'], nPairCorr);
    fprintf(fid, '\\midrule\n\n');

    condLabelMap = containers.Map( ...
        {'original','gatys','texform'}, ...
        {'Natural','Gatys','Texform'});

    for r = 1:nROI
        roi_name = rois{r};
        nRowsROI = nCond;   % one row per condition

        for c = 1:nCond
            cond_name = conditions{c};
            rowMask = strcmp(T_pair_all.ROI, roi_name) & strcmp(T_pair_all.Condition, cond_name);
            row = T_pair_all(rowMask, :);

            if isempty(row)
                continue;
            end

            if c == 1
                fprintf(fid, '\\multirow{%d}{*}{\\textbf{%s}} ', nRowsROI, localPrettyROI(roi_name));
            else
                fprintf(fid, ' ');
            end

            if row.TargetIsHighest
                targetFlag = 'yes';
            else
                targetFlag = 'no';
            end

            fprintf(fid, '& %s & %s & %s & %s & %d & %.2f & %d & %s & %s \\\\\n', ...
                condLabelMap(cond_name), ...
                localPrettyCat(char(row.Cat1)), ...
                localPrettyCat(char(row.Cat2)), ...
                targetFlag, ...
                row.n, ...
                row.t, ...
                row.df, ...
                format_p_label_plain(row.p_uncorr), ...
                format_p_label_plain(row.p_bonf));
        end

        if r < nROI
            fprintf(fid, '\n\\midrule\n\n');
        end
    end

    fprintf(fid, '\n\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');

    fclose(fid);
    fprintf('Saved LaTeX pairwise table: %s\n', latexFilePair);
end

%% ============================================================
%% LATEX EXPORT: BASELINE (VS ZERO) SIGNIFICANCE TABLE
%% ============================================================
% Rows: ROI x Condition x Category
% Bonferroni N = nBaseCorr = 15 (3 ROIs x 5 categories, applied per condition)
if ~isempty(T_baseline_all)

    latexFileBase = fullfile(outputDir, 'fmri_all_rois_baseline_table.tex');
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
        {'original','gatys','texform'}, ...
        {'Natural','Gatys','Texform'});

    for r = 1:nROI
        roi_name = rois{r};
        nRowsROI = nCond * nCat;

        first_roi_row = true;
        for c = 1:nCond
            cond_name = conditions{c};
            nRowsCond = nCat;

            first_cond_row = true;
            for g = 1:nCat
                this_cat = categories{g};
                rowMask = strcmp(T_baseline_all.ROI, roi_name) ...
                        & strcmp(T_baseline_all.Condition, cond_name) ...
                        & strcmp(T_baseline_all.Category, this_cat);
                row = T_baseline_all(rowMask, :);

                if isempty(row)
                    continue;
                end

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
                    localPrettyCat(this_cat), ...
                    row.n, ...
                    row.mean, ...
                    row.sd, ...
                    row.t, ...
                    row.df, ...
                    format_p_label_plain(row.p_uncorr), ...
                    format_p_label_plain(row.p_bonf));
            end

            % thin rule between conditions (except after last)
            if c < nCond
                fprintf(fid, '\\cmidrule(l){2-10}\n');
            end
        end

        if r < nROI
            fprintf(fid, '\n\\midrule\n\n');
        end
    end

    fprintf(fid, '\n\\bottomrule\n');
    fprintf(fid, '\\end{tabular}}\n');
    fprintf(fid, '\\end{table}\n');

    fclose(fid);
    fprintf('Saved LaTeX baseline table: %s\n', latexFileBase);
end

%% ============================================================
%% HELPER FUNCTIONS
%% ============================================================
function test_vs0 = test_against_zero(metric_vals, n_cat, bonf_n)

    p_vals  = nan(1, n_cat);
    t_vals  = nan(1, n_cat);
    df_vals = nan(1, n_cat);

    for i = 1:n_cat
        vals = metric_vals{i};
        vals = vals(~isnan(vals));

        if numel(vals) < 2
            continue;
        end

        [~, p_vals(i), ~, stats] = ttest(vals, 0);
        t_vals(i)  = stats.tstat;
        df_vals(i) = stats.df;
    end

    test_vs0.p = p_vals;
    test_vs0.p_bonf = min(p_vals * bonf_n, 1);   % corrected for bonf_n tests
    test_vs0.t = t_vals;
    test_vs0.df = df_vals;
end

function pair_test = test_highest_vs_second(metric_vals, metric_mean, categories, targetCatIdx)
% Compare the target ROI category against the highest-d' category.
% If the target IS already the highest, fall back to highest vs second-highest.


    [~, ord] = sort(metric_mean, 'descend', 'MissingPlacement', 'last');

    if numel(ord) < 2 || any(isnan(metric_mean(ord(1:2))))
        pair_test.idx1            = NaN;
        pair_test.idx2            = NaN;
        pair_test.cat1            = '';
        pair_test.cat2            = '';
        pair_test.isTargetHighest = false;
        pair_test.n               = NaN;
        pair_test.t               = NaN;
        pair_test.df              = NaN;
        pair_test.p               = NaN;
        return;
    end

    highestIdx = ord(1);

    if highestIdx == targetCatIdx
        % Target is already the highest: compare highest vs second-highest
        idx1 = highestIdx;
        idx2 = ord(2);
        isTargetHighest = true;
    else
        % Compare highest vs target
        idx1 = highestIdx;
        idx2 = targetCatIdx;
        isTargetHighest = false;
    end

    v1 = metric_vals{idx1};
    v2 = metric_vals{idx2};

    validMask = ~isnan(v1) & ~isnan(v2);
    v1 = v1(validMask);
    v2 = v2(validMask);

    if numel(v1) < 2
        tstat  = NaN;
        df     = NaN;
        p_p    = NaN;
        n_pair = numel(v1);
    else
        [~, p_p, ~, stats] = ttest(v1, v2);
        tstat  = stats.tstat;
        df     = stats.df;
        n_pair = numel(v1);
    end

    pair_test.idx1            = idx1;
    pair_test.idx2            = idx2;
    pair_test.cat1            = categories{idx1};
    pair_test.cat2            = categories{idx2};
    pair_test.isTargetHighest = isTargetHighest;
    pair_test.n               = n_pair;
    pair_test.t               = tstat;
    pair_test.df              = df;
    pair_test.p               = p_p;
end

function save_metric_csvs(outDir, file_tag, categories, n_vals, m_mean, m_sem, m_sd, test_vs0, pair_test)

    % --- per-category summary stats ---
    T_stats = table( ...
        categories(:), ...
        n_vals(:), ...
        m_mean(:), ...
        m_sem(:), ...
        m_sd(:), ...
        test_vs0.t(:), ...
        test_vs0.df(:), ...
        test_vs0.p(:), ...
        test_vs0.p_bonf(:), ...
        'VariableNames', { ...
            'category', ...
            'n', ...
            'mean', ...
            'sem', ...
            'sd', ...
            't_vs_0', ...
            'df_vs_0', ...
            'p_vs_0', ...
            'p_bonf_vs_0'} ...
    );

    writetable(T_stats, fullfile(outDir, [file_tag '_stats.csv']));

    % --- highest vs target pairwise test ---
    T_pair = table( ...
        string(pair_test.cat1), ...
        string(pair_test.cat2), ...
        pair_test.isTargetHighest, ...
        pair_test.idx1, ...
        pair_test.idx2, ...
        pair_test.n, ...
        pair_test.t, ...
        pair_test.df, ...
        pair_test.p, ...
        pair_test.p_bonf, ...
        'VariableNames', { ...
            'cat1_highest', ...
            'cat2_target_or_second', ...
            'target_is_highest', ...
            'idx_cat1', ...
            'idx_cat2', ...
            'n_paired', ...
            't_cat1_vs_cat2', ...
            'df_cat1_vs_cat2', ...
            'p_cat1_vs_cat2', ...
            'p_bonf_cat1_vs_cat2'} ...
    );

    writetable(T_pair, fullfile(outDir, [file_tag '_pairwise_target_vs_highest.csv']));
end

function plot_metric_barplot(m_mean, m_sem, p_bonf, pair_test, categories, bar_colors, clean_name, ylab, outFile, yLims)
    n_cat = numel(categories);

    fig = figure('Position', [100 100 700 520], 'Visible', 'off');
    hold on;

    b = bar(1:n_cat, m_mean, 'FaceColor', 'flat', 'EdgeColor', 'none');
    for c = 1:n_cat
        b.CData(c,:) = bar_colors(c,:);
    end

    errorbar(1:n_cat, m_mean, m_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

    ylim(yLims);
    y_span = range(ylim);
    star_offset = 0.05 * y_span;

    % Per-category baseline significance
   
    for i = 1:n_cat
        if isnan(p_bonf(i)) || isnan(m_mean(i)) || isnan(m_sem(i))
            continue; 
        end

        if m_mean(i) < 0
            continue;  
        end

        if p_bonf(i) < 0.001
            s = '***';
        elseif p_bonf(i) < 0.01
            s = '**';
        elseif p_bonf(i) < 0.05
            s = '*';
        else
            s = 'n.s.';
        end

        y_text = m_mean(i) + m_sem(i) + star_offset;
        text(i, y_text, s, 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    end

    % Pairwise bracket significance (uses Bonferroni-corrected p)
    if isfield(pair_test, 'p_bonf') && ...
            ~isnan(pair_test.p_bonf) && ~isnan(pair_test.idx1) && ~isnan(pair_test.idx2) && ...
            ~isnan(m_mean(pair_test.idx1)) && ~isnan(m_mean(pair_test.idx2))

        p_pair_corr = pair_test.p_bonf;

        if p_pair_corr < 0.001
            star_pair = '***';
        elseif p_pair_corr < 0.01
            star_pair = '**';
        elseif p_pair_corr < 0.05
            star_pair = '*';
        else
            star_pair = 'n.s.';
        end

        y_base = max(m_mean + m_sem, [], 'omitnan');
        pad = 0.08 * y_span;
        h = 0.03 * y_span;

        plot([pair_test.idx1 pair_test.idx1 pair_test.idx2 pair_test.idx2], ...
             [y_base+pad y_base+pad+h y_base+pad+h y_base+pad], ...
             'k-', 'LineWidth', 1.5);

        text(mean([pair_test.idx1 pair_test.idx2]), y_base+pad+h+0.03*y_span, star_pair, ...
            'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
    end

    ylabel(ylab);
    title(clean_name, 'Interpreter', 'none');
    box off;
    set(gca, 'XTick', 1:n_cat, 'XTickLabel', categories);

    exportgraphics(fig, outFile, 'Resolution', 300);
    %close(fig);
end

function pretty_roi = localPrettyROI(roi_name)
    switch roi_name
        case 'LOTC_body'
            pretty_roi = 'LOTC Body';
        case 'LOTC_hand'
            pretty_roi = 'LOTC Hand';
        case 'LOTC_tool'
            pretty_roi = 'LOTC Tool';
        otherwise
            pretty_roi = strrep(roi_name, '_', ' ');
    end
end

function pretty_cat = localPrettyCat(cat_name)
    switch cat_name
        case 'body'
            pretty_cat = 'Bodies';
        case 'hand'
            pretty_cat = 'Hands';
        case 'tool'
            pretty_cat = 'Tools';
        case 'manip'
            pretty_cat = 'Manip';
        case 'nonmanip'
            pretty_cat = 'Nman';
        otherwise
            pretty_cat = cat_name;
    end
end

function p_str = format_p_label(p)
  
    if isnan(p)
        p_str = 'n/a';
    elseif p < 0.001
        p_str = '$p {<} .001$';
    else
        % round values to 3 decimals
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