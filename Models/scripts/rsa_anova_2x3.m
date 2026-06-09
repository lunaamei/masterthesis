
% 2 (Model: baseline, dvdmodel) x 3 (Condition: natural, gatys, texform)
% repeated-measures ANOVA + one-sample t-tests + condition comparisons
% + model comparisons, per ROI.
%
% INPUT
%   per_subject_rho.csv
%
%
% OUTPUTS 
%   all_rois_rmanova.csv    
%   all_rois_onesample.csv  
%   all_rois_cond_within_model.csv 
%   all_rois_model_comparison.csv 
%   all_rois_descriptives.csv     
%   rsa_anova_table.tex         
%   rsa_comprehensive_table.tex   
%   rsa_cond_within_model_table.tex

clear; clc;

%% ============================================================
%% PATHS & SETTINGS
%% ============================================================
inputCsv  = '/Users/lunameidoering/Desktop/CODE/Models/out/rdm_13subs/per_subject_rho.csv';

outputDir = '/Users/lunameidoering/Desktop/CODE/Models/out/rdm_13subs/anova/';

if ~exist(outputDir, 'dir'), mkdir(outputDir); end

conditions = {'natural','gatys','texform'};
models     = {'baseline','dvdmodel'};
rois       = {'body','hand','tool'};
roiDisp    = {'LOTC Body Left','LOTC Hand Left','LOTC Tool Left'};
condDisp   = {'Natural','Gatys','Texform'};
modelDisp  = {'Baseline','DVD Model'};

%% ============================================================
%% LOAD & VALIDATE
%% ============================================================
if ~exist(inputCsv,'file')
    error('per_subject_rho.csv not found:\n  %s\nRun the Python notebook first.', inputCsv);
end

T = readtable(inputCsv, 'Delimiter', ';');
for col = {'subject','model','condition','roi'}
    if ~ismember(col{1}, T.Properties.VariableNames)
        error('Missing column: %s', col{1});
    end
    if ~iscellstr(T.(col{1}))
        T.(col{1}) = cellstr(T.(col{1}));
    end
end
fprintf('Loaded %d rows from per_subject_rho.csv\n', height(T));

%% ============================================================
%% PER-ROI ANALYSIS
%% ============================================================


all_anova_rows  = {};  
all_os_rows     = {};   
all_cond_rows   = {};   
all_comp_rows   = {};   
all_desc_rows   = {};  

wideMatAll = cell(numel(rois), 1);
withinAll  = cell(numel(rois), 1);
nSubjAll   = zeros(numel(rois), 1);

for r = 1:numel(rois)
    thisRoi = rois{r};
    Troi     = T(strcmp(T.roi, thisRoi), :);
    subjects = unique(Troi.subject, 'stable');
    nSubjIn  = numel(subjects);

    % Build wide matrix
    nMeas      = numel(models) * numel(conditions);  % 6
    measNames  = cell(1, nMeas);
    withinMod  = cell(1, nMeas);
    withinCond = cell(1, nMeas);
    wideMat    = nan(nSubjIn, nMeas);

    k = 0;
    for mi = 1:numel(models)
        for ci = 1:numel(conditions)
            k = k + 1;
            measNames{k}   = sprintf('%s_%s', models{mi}, conditions{ci});
            withinMod{k}   = models{mi};
            withinCond{k}  = conditions{ci};
            for si = 1:nSubjIn
                mask = strcmp(Troi.subject, subjects{si}) & ...
                       strcmp(Troi.model,   models{mi})   & ...
                       strcmp(Troi.condition, conditions{ci});
                if any(mask)
                    wideMat(si, k) = Troi.rho(find(mask, 1));
                end
            end
        end
    end

    % Drop subjects with any missing cell 
    completeMask = all(~isnan(wideMat), 2);
    if any(~completeMask)
        fprintf('  Dropping %d subject(s): %s\n', sum(~completeMask), ...
            strjoin(subjects(~completeMask), ', '));
    end
    wideMat  = wideMat(completeMask, :);
    subjects = subjects(completeMask);
    nSubj    = numel(subjects);

    wideMatAll{r} = wideMat;
    withinAll{r}  = struct('model', {withinMod}, 'cond', {withinCond});
    nSubjAll(r)   = nSubj;

    % Descriptives
    for k = 1:nMeas
        v = wideMat(:, k);
        all_desc_rows(end+1, :) = {thisRoi, withinMod{k}, withinCond{k}, ...
            nSubj, mean(v), std(v), std(v)/sqrt(nSubj)};
    end

    % ANOVA 
    wideTbl = array2table(wideMat, 'VariableNames', measNames);
    wideTbl = addvars(wideTbl, subjects, 'Before', 1, 'NewVariableNames', 'Subject');

    withinDesign = table( ...
        categorical(withinMod',  models,     'Ordinal', true), ...
        categorical(withinCond', conditions, 'Ordinal', true), ...
        'VariableNames', {'Model','Condition'});

    measRange = sprintf('%s-%s ~ 1', measNames{1}, measNames{end});
    rm  = fitrm(wideTbl, measRange, 'WithinDesign', withinDesign);
    aov = ranova(rm, 'WithinModel', 'Model*Condition');
    disp(aov);
    writetable(aov, fullfile(outputDir, sprintf('ranova_%s.csv', thisRoi)), ...
        'WriteRowNames', true, 'Delimiter', ';');

    rn = aov.Properties.RowNames;
    for effCell = {'Model','Condition','Model:Condition'}
        eff     = effCell{1};
        rowName = ['(Intercept):' eff];
        errName = ['Error(' eff ')'];
        if ~any(strcmp(rn, rowName)), continue; end
        F     = aov.F(strcmp(rn, rowName));
        df1   = aov.DF(strcmp(rn, rowName));
        df2   = aov.DF(strcmp(rn, errName));
        p     = aov.pValue(strcmp(rn, rowName));
        pGG   = aov.pValueGG(strcmp(rn, rowName));
        ssE   = aov.SumSq(strcmp(rn, rowName));
        ssR   = aov.SumSq(strcmp(rn, errName));
        pEta2 = ssE / (ssE + ssR);
        all_anova_rows(end+1, :) = {thisRoi, eff, F, df1, df2, p, pGG, pEta2, nSubj};
    end

    % one-sample t-tests 
    % 2 models × 3 conditions = 6 tests per ROI.
    n_os      = numel(models) * numel(conditions);  % 6
    os_p_raw  = zeros(1, n_os);
    os_t      = zeros(1, n_os);
    os_df_    = zeros(1, n_os);
    os_mean   = zeros(1, n_os);
    os_mod    = cell(1, n_os);
    os_cond_  = cell(1, n_os);

    k = 0;
    for mi = 1:numel(models)
        for ci = 1:numel(conditions)
            k = k + 1;
            col  = strcmp(withinMod, models{mi}) & strcmp(withinCond, conditions{ci});
            vals = wideMat(:, col);
            [~, os_p_raw(k), ~, st] = ttest(vals, 0);
            os_t(k)    = st.tstat;
            os_df_(k)  = st.df;
            os_mean(k) = mean(vals);
            os_mod{k}  = models{mi};
            os_cond_{k} = conditions{ci};
        end
    end

    os_p_bonf = min(os_p_raw * n_os, 1); 

    
    for k = 1:n_os
        fprintf('  %-10s %-8s: t(%d)=%6.3f,  p_raw=%.4g,  p_bonf(%d)=%.4g\n', ...
            os_mod{k}, os_cond_{k}, os_df_(k), os_t(k), ...
            os_p_raw(k), n_os, os_p_bonf(k));
        all_os_rows(end+1, :) = {thisRoi, os_mod{k}, os_cond_{k}, nSubj, ...
            os_mean(k), os_t(k), os_df_(k), os_p_raw(k), os_p_bonf(k)};
    end

 %within-model condition comparisons 
% 2 models × 3 pairs = 6 tests per ROI, all corrected together.

pairs        = nchoosek(1:numel(conditions), 2); 
n_fam2       = numel(models) * size(pairs, 1); 

% Pre-allocate collectors for all 6 tests
f2_p_raw  = zeros(1, n_fam2);
f2_t      = zeros(1, n_fam2);
f2_df     = zeros(1, n_fam2);
f2_mod    = cell(1, n_fam2);
f2_c1     = cell(1, n_fam2);
f2_c2     = cell(1, n_fam2);

idx = 0;
for mi = 1:numel(models)
    mod_mat = nan(nSubj, numel(conditions));
    for ci = 1:numel(conditions)
        col = strcmp(withinMod, models{mi}) & strcmp(withinCond, conditions{ci});
        mod_mat(:, ci) = wideMat(:, col);
    end
    for pp = 1:size(pairs, 1)       
        idx = idx + 1;
        [~, f2_p_raw(idx), ~, st] = ttest(mod_mat(:, pairs(pp,1)), ...
                                           mod_mat(:, pairs(pp,2)));
        f2_t(idx)    = st.tstat;
        f2_df(idx)   = st.df;
        f2_mod{idx}  = models{mi};
        f2_c1{idx}   = conditions{pairs(pp,1)};
        f2_c2{idx}   = conditions{pairs(pp,2)};
    end
end

% Apply Bonferroni × 6 to all 6 p-values together
f2_p_bonf = min(f2_p_raw * n_fam2, 1);

fprintf('\n--- Family 2: within-model condition comparisons [Bonferroni x%d] ---\n', n_fam2);
for k = 1:n_fam2
    fprintf('  %-10s  %s vs %-8s: t(%d)=%6.3f,  p_raw=%.4g,  p_bonf(%d)=%.4g\n', ...
        f2_mod{k}, f2_c1{k}, f2_c2{k}, f2_df(k), f2_t(k), ...
        f2_p_raw(k), n_fam2, f2_p_bonf(k));
    all_cond_rows(end+1, :) = {thisRoi, f2_mod{k}, f2_c1{k}, f2_c2{k}, nSubj, ...
        f2_t(k), f2_df(k), f2_p_raw(k), f2_p_bonf(k)};
end

    % model comparison within condition 
    nConds    = numel(conditions);  % 3
    mc_p_raw  = zeros(1, nConds);
    mc_t_     = zeros(1, nConds);
    mc_df_    = zeros(1, nConds);

    for ci = 1:numel(conditions)
        c_base = strcmp(withinCond, conditions{ci}) & strcmp(withinMod, 'baseline');
        c_dvd  = strcmp(withinCond, conditions{ci}) & strcmp(withinMod, 'dvdmodel');
        [~, mc_p_raw(ci), ~, st] = ttest(wideMat(:,c_dvd), wideMat(:,c_base));
        mc_t_(ci)  = st.tstat;
        mc_df_(ci) = st.df;
    end

    mc_p_bonf = min(mc_p_raw * nConds, 1);   % Bonferroni x3

    fprintf('\n--- Family 3: DVD - Baseline within condition [Bonferroni x%d] ---\n', nConds);
    for ci = 1:numel(conditions)
        fprintf('  %-8s: t(%d)=%6.3f,  p_raw=%.4g,  p_bonf(%d)=%.4g\n', ...
            conditions{ci}, mc_df_(ci), mc_t_(ci), ...
            mc_p_raw(ci), nConds, mc_p_bonf(ci));
        all_comp_rows(end+1, :) = {thisRoi, conditions{ci}, nSubj, ...
            mc_t_(ci), mc_df_(ci), mc_p_raw(ci), mc_p_bonf(ci)};
    end
end

%% ============================================================
%% 4. SAVE TABLES
%% ============================================================
T_anova = cell2table(all_anova_rows, 'VariableNames', ...
    {'ROI','Effect','F','df1','df2','p','pGG','partial_eta2','n'});

T_onesample = cell2table(all_os_rows, 'VariableNames', ...
    {'ROI','Model','Condition','n','mean_rho','t','df','p_raw','p_bonf_6'});

T_cond_mod = cell2table(all_cond_rows, 'VariableNames', ...
    {'ROI','Model','Cond1','Cond2','n','t','df','p_raw','p_bonf_6'});

T_model_comp = cell2table(all_comp_rows, 'VariableNames', ...
    {'ROI','Condition','n','t','df','p_raw','p_bonf_3'});

T_desc = cell2table(all_desc_rows, 'VariableNames', ...
    {'ROI','Model','Condition','n','mean','sd','sem'});

writetable(T_anova,      fullfile(outputDir,'all_rois_rmanova.csv'),            'Delimiter',';');
writetable(T_onesample,  fullfile(outputDir,'all_rois_onesample.csv'),          'Delimiter',';');
writetable(T_cond_mod,   fullfile(outputDir,'all_rois_cond_within_model.csv'),  'Delimiter',';');
writetable(T_model_comp, fullfile(outputDir,'all_rois_model_comparison.csv'),   'Delimiter',';');
writetable(T_desc,       fullfile(outputDir,'all_rois_descriptives.csv'),       'Delimiter',';');

fprintf('\n=== ANOVA ===\n');        disp(T_anova);
fprintf('=== Family 1 (one-sample) ===\n'); disp(T_onesample);
fprintf('=== Family 2 (cond within model) ===\n'); disp(T_cond_mod);
fprintf('=== Family 3 (model comparison) ===\n');  disp(T_model_comp);

%% ============================================================
%% 6. LATEX EXPORT
%% ============================================================
% Three tables:
%   Table 1: rmANOVA effects.
%   Table 2: Comprehensive results – one-sample (Family 1) and model
%            comparison (Family 3), both p_raw and p_bonf shown.
%   Table 3: Within-model condition comparisons (Family 2).

effect_display = containers.Map( ...
    {'Model','Condition','Model:Condition'}, ...
    {'Model','Condition','Model $\times$ Condition'});

%% ---- Table 1: rmANOVA ─────────────────────────────────────────────────────
latexAnova = fullfile(outputDir, 'rsa_anova_table.tex');
fid = fopen(latexAnova, 'w');
fprintf(fid, '%% Auto-generated: rsa_anova_from_python.m\n');
fprintf(fid, '\\begin{table}[!ht]\n\\centering\n');
fprintf(fid, '\\begin{tabular}{llccccc}\n\\toprule\n\n');
fprintf(fid, '\\textbf{ROI} & \\textbf{Effect} & \\textbf{F} & \\textbf{df1} & \\textbf{df2} & \\textbf{$p_{\\text{GG}}$} & \\textbf{$\\eta_p^2$} \\\\\n\n');
fprintf(fid, '\\midrule\n\n');

eff_order = {'Model','Condition','Model:Condition'};
for r = 1:numel(rois)
    roi_name  = rois{r};
    roi_label = roiDisp{strcmp(rois, roi_name)};
    T_roi = T_anova(strcmp(T_anova.ROI, roi_name), :);
    [~, ord] = ismember(T_roi.Effect, eff_order);
    [~, srt] = sort(ord);  T_roi = T_roi(srt, :);
    for i = 1:height(T_roi)
        if i == 1
            fprintf(fid, '\\multirow{3}{*}{\\textbf{%s}} ', roi_label);
        else
            fprintf(fid, ' ');
        end
        eff_lbl = effect_display(T_roi.Effect{i});
        fprintf(fid, '& %s & %.2f & %d & %d & %s & %.3f \\\\\n', ...
            eff_lbl, T_roi.F(i), T_roi.df1(i), T_roi.df2(i), ...
            fmt_p_exact(T_roi.pGG(i)), T_roi.partial_eta2(i));
    end
    if r < numel(rois), fprintf(fid, '\n\\midrule\n\n'); end
end
fprintf(fid, '\n\\bottomrule\n\\end{tabular}\n');
fprintf(fid, '\\caption{Repeated-measures ANOVA (2 Models $\\times$ 3 Conditions) per ROI.\n');
fprintf(fid, '$p_{\\text{GG}}$~= Greenhouse--Geisser corrected $p$.\n');
fprintf(fid, '$\\eta_p^2$~= partial eta squared.\n');
fprintf(fid, '$^{*}~p<.05$, $^{**}~p<.01$, $^{***}~p<.001$.}\n');
fprintf(fid, '\\label{tab:anova}\n\\end{table}\n');
fclose(fid);
fprintf('Saved: %s\n', latexAnova);

%% ---- Table 2: Comprehensive ─────────────────────────

latexComp = fullfile(outputDir, 'rsa_comprehensive_table.tex');
fid = fopen(latexComp, 'w');
fprintf(fid, '%% Auto-generated: rsa_anova_from_python.m\n');
fprintf(fid, '\\begin{table}[!ht]\n\\centering\n');
fprintf(fid, '\\resizebox{\\textwidth}{!}{%%\n');
fprintf(fid, '\\begin{tabular}{ll cccc @{\\hspace{0.6em}} cccc @{\\hspace{0.6em}} ccc}\n');
fprintf(fid, '\\toprule\n\n');

fprintf(fid, '\\multirow{2}{*}{\\textbf{ROI}} & \\multirow{2}{*}{\\textbf{Condition}} ');
fprintf(fid, '& \\multicolumn{4}{c}{\\textbf{DVD-shape}} ');
fprintf(fid, '& \\multicolumn{4}{c}{\\textbf{DVD-texture}} ');
fprintf(fid, '& \\multicolumn{3}{c}{\\textbf{Model Difference}} \\\\\n\n');
fprintf(fid, '\\cmidrule(lr){3-6} \\cmidrule(lr){7-10} \\cmidrule(lr){11-13}\n\n');
fprintf(fid, ' & & $r$ & $t(df)$ & $p$ & $p_{\\text{bonf}(6)}$ ');
fprintf(fid, '& $r$ & $t(df)$ & $p$ & $p_{\\text{bonf}(6)}$ ');
fprintf(fid, '& $t(df)$ & $p$ & $p_{\\text{bonf}(3)}$ \\\\\n\n');
fprintf(fid, '\\midrule\n\n');

for r = 1:numel(rois)
    roi_name  = rois{r};
    roi_label = roiDisp{strcmp(rois, roi_name)};

    for ci = 1:numel(conditions)
        cond_name  = conditions{ci};
        cond_label = condDisp{strcmp(conditions, cond_name)};

        row_dvd  = T_onesample(strcmp(T_onesample.ROI, roi_name) & ...
                               strcmp(T_onesample.Model, 'dvdmodel') & ...
                               strcmp(T_onesample.Condition, cond_name), :);
        row_base = T_onesample(strcmp(T_onesample.ROI, roi_name) & ...
                               strcmp(T_onesample.Model, 'baseline') & ...
                               strcmp(T_onesample.Condition, cond_name), :);
        row_comp = T_model_comp(strcmp(T_model_comp.ROI, roi_name) & ...
                                strcmp(T_model_comp.Condition, cond_name), :);

        dvd_r  = fmt_r_val(row_dvd.mean_rho);
        dvd_t  = fmt_t_val(row_dvd.t, row_dvd.df);
        dvd_p  = fmt_p_exact(row_dvd.p_raw);
        dvd_pb = fmt_p_exact(row_dvd.p_bonf_6);

        base_r  = fmt_r_val(row_base.mean_rho);
        base_t  = fmt_t_val(row_base.t, row_base.df);
        base_p  = fmt_p_exact(row_base.p_raw);
        base_pb = fmt_p_exact(row_base.p_bonf_6);

        comp_t  = fmt_t_val(row_comp.t, row_comp.df);
        comp_p  = fmt_p_exact(row_comp.p_raw);
        comp_pb = fmt_p_exact(row_comp.p_bonf_3);

        if ci == 1
            roi_str = sprintf('\\multirow{3}{*}{\\textbf{%s}}', roi_label);
        else
            roi_str = ' ';
        end

        fprintf(fid, '%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\\n', ...
            roi_str, cond_label, ...
            dvd_r,  dvd_t,  dvd_p,  dvd_pb, ...
            base_r, base_t, base_p, base_pb, ...
            comp_t, comp_p, comp_pb);
    end
    if r < numel(rois), fprintf(fid, '\n\\midrule\n\n'); end
end

fprintf(fid, '\n\\bottomrule\n\\end{tabular}}\n');
fprintf(fid, '\\caption{One-sample $t$-tests against zero and model comparison per ROI\n');
fprintf(fid, 'and condition. $r$~= mean Pearson correlation (model--fMRI). $t(df)$~=\n');
fprintf(fid, '$t$-statistic with degrees of freedom. $p$~= uncorrected $p$-value.\n');
fprintf(fid, '$p_{\\text{bonf}(6)}$~= Bonferroni-corrected $p$ for the one-sample family\n');
fprintf(fid, '(6 tests per ROI); $p_{\\text{bonf}(3)}$~= Bonferroni-corrected $p$ for the\n');
fprintf(fid, 'model comparison family (3 tests per ROI).\n');
fprintf(fid, '$^{*}~p<.05$, $^{**}~p<.01$, $^{***}~p<.001$ (based on uncorrected $p$).}\n');
fprintf(fid, '\\label{tab:comprehensive}\n\\end{table}\n');
fclose(fid);
fprintf('Saved: %s\n', latexComp);

%% ---- Table 3: Within-model condition comparisons (Family 2) ────────────────
latexCond = fullfile(outputDir, 'rsa_cond_within_model_table.tex');
fid = fopen(latexCond, 'w');
fprintf(fid, '%% Auto-generated: rsa_anova_from_python.m\n');
fprintf(fid, '\\begin{table}[!ht]\n\\centering\n');
fprintf(fid, '\\begin{tabular}{lll cccc}\n\\toprule\n\n');
fprintf(fid, '\\textbf{ROI} & \\textbf{Model} & \\textbf{Comparison} & ');
fprintf(fid, '\\textbf{$t(df)$} & \\textbf{$p$} & \\textbf{$p_{\\text{bonf}(3)}$} \\\\\n\n');
fprintf(fid, '\\midrule\n\n');

pair_labels = containers.Map( ...
    {'natural_gatys','natural_texform','gatys_texform'}, ...
    {'Natural -- Gatys','Natural -- Texform','Gatys -- Texform'});

pairs = nchoosek(1:numel(conditions), 2);

for r = 1:numel(rois)
    roi_name  = rois{r};
    roi_label = roiDisp{strcmp(rois, roi_name)};
    n_rows_roi = numel(models) * size(pairs,1);  % 6
    row_idx = 0;

    for mi = 1:numel(models)
        mod_label = modelDisp{strcmp(models, models{mi})};
        n_rows_mod = size(pairs, 1);  % 3

        for pp = 1:size(pairs,1)
            row_idx = row_idx + 1;
            c1 = conditions{pairs(pp,1)};
            c2 = conditions{pairs(pp,2)};
            pair_key   = sprintf('%s_%s', c1, c2);
            pair_label = pair_labels(pair_key);

            row = T_cond_mod(strcmp(T_cond_mod.ROI,   roi_name) & ...
                             strcmp(T_cond_mod.Model,  models{mi}) & ...
                             strcmp(T_cond_mod.Cond1,  c1) & ...
                             strcmp(T_cond_mod.Cond2,  c2), :);

            if row_idx == 1
                roi_str = sprintf('\\multirow{%d}{*}{\\textbf{%s}}', n_rows_roi, roi_label);
            else
                roi_str = ' ';
            end
            if pp == 1
                mod_str = sprintf('\\multirow{%d}{*}{%s}', n_rows_mod, mod_label);
            else
                mod_str = ' ';
            end

            fprintf(fid, '%s & %s & %s & %s & %s & %s \\\\\n', ...
                roi_str, mod_str, pair_label, ...
                fmt_t_val(row.t, row.df), ...
                fmt_p_exact(row.p_raw), ...
                fmt_p_exact(row.p_bonf_6));
        end
        if mi < numel(models)
            fprintf(fid, '\\cmidrule(lr){2-6}\n');
        end
    end
    if r < numel(rois), fprintf(fid, '\n\\midrule\n\n'); end
end

fprintf(fid, '\n\\bottomrule\n\\end{tabular}\n');
fprintf(fid, '\\caption{Within-model pairwise condition comparisons (Family~2).\n');
fprintf(fid, 'Paired $t$-tests; Bonferroni-corrected across 3 pairs per model per ROI.\n');
fprintf(fid, '$p$~= uncorrected; $p_{\\text{bonf}(3)}$~= Bonferroni-corrected.\n');
fprintf(fid, '$^{*}~p<.05$, $^{**}~p<.01$, $^{***}~p<.001$ (based on uncorrected $p$).}\n');
fprintf(fid, '\\label{tab:cond_within_model}\n\\end{table}\n');
fclose(fid);
fprintf('Saved: %s\n', latexCond);

fprintf('\nDone. Outputs in: %s\n', outputDir);

%% ============================================================
%% HELPERS
%% ============================================================
function s = fmt_p_exact(p_raw)
    % Exact p to 3 decimal places with significance stars (* <.05 ** <.01 *** <.001).
    if isnan(p_raw),   s = '$-$';     return; end
    if p_raw >= 1,     s = '$1.000$'; return; end
    if p_raw < 0.001
        val = '<.001';  stars = '^{***}';
    else
        tmp  = sprintf('%.3f', p_raw);
        val  = tmp(2:end);   % strip leading zero
        if     p_raw < 0.01, stars = '^{**}';
        elseif p_raw < 0.05, stars = '^{*}';
        else,                stars = '';
        end
    end
    s = sprintf('$%s%s$', val, stars);
end

function s = fmt_r_val(r_val)
    s = sprintf('$%.2f$', r_val);
    s = strrep(s, '$0.', '$.');
    s = strrep(s, '$-0.', '$-.');
end

function s = fmt_t_val(t_val, df_val)
    s = sprintf('$%.2f$ (%.0f)', t_val, df_val);
end
