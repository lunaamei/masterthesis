
% updated for using CV!
% THIS IS THE LATEST VERSION!!!!
%
% For each prediction file:
% - loads cv_predictions and all_fullpaths_sorted
% - assigns stimuli to categories
% - computes raw activation statistics
% - computes per-stimulus d' values using pooled variance:
%       d'_stim = (x_cat - mean(noncat)) / sqrt(0.5*(var(cat)+var(noncat)))
% - saves:
%       raw_activation_stats.csv
%       raw_activation_pairwise_highest_vs_second.csv
%       dprime_selectivity_stats.csv
%       dprime_selectivity_pairwise_highest_vs_second.csv
%       dprime_selectivity_single_stimulus.csv
%       raw_activation_barplot.png
%       dprime_selectivity_barplot.png
% =========================================================================

clear; clc; close all;

%% ===== PATHS =====

% Baseline ImageNet
%predDir = '/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/EncodingModels/final/baseline_imagenet_cv/predictions_cv/';
%baseOutDir = '/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/EncodingModels/final/baseline_imagenet_cv/results/';


% DVD ImageNet
predDir = '/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/EncodingModels/final/dvd_imagenet_cv/predictions_cv/';
baseOutDir = '/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/EncodingModels/final/dvd_imagenet_cv/results/';

if ~exist(baseOutDir, 'dir')
    mkdir(baseOutDir);
end

addpath(predDir);

%% ===== FIND ALL .MAT FILES =====
predFiles = dir(fullfile(predDir, '*.mat'));

if isempty(predFiles)
    error('No .mat files found in %s', predDir);
end

fprintf('Found %d prediction files.\n', numel(predFiles));

%% ===== CATEGORY DEFINITIONS =====
bodies_stim = [1:12,   101:113];
hands_stim  = [13:25,  114:125];
tools_stim  = [26:50,  126:150];
mani_stim   = [51:75,  151:175];
nman_stim   = [76:100, 176:200];

categories    = {'bodies','hands','tools','mani','nman'};
category_stim = {bodies_stim, hands_stim, tools_stim, mani_stim, nman_stim};
n_cat         = numel(categories);

%% ===== COLORS =====
bar_colors = [
    1, 0.5, 0;        % bodies
    1, 0.9, 0;        % hands
    0, 0, 0.5;        % tools
    0.68, 0.85, 0.9;  % mani
    0.75, 0, 0.75     % nman
];

%% =================================================================
%                      LOOP OVER ALL FILES
%% =================================================================
for f = 1:numel(predFiles)

    pred_name = predFiles(f).name;
    pred_path = fullfile(predDir, pred_name);

    fprintf('\n============================================================\n');
    fprintf('Processing file %d / %d\n', f, numel(predFiles));
    fprintf('%s\n', pred_name);
    fprintf('============================================================\n');

    [~, clean_name, ~] = fileparts(pred_name);
    clean_name_sp = strrep(clean_name, '_', ' ');

    outDir = fullfile(baseOutDir, pred_name);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    try
        %% ===== LOAD PREDICTIONS =====
        predictions = load(pred_path);

        if ~isfield(predictions, 'cv_predictions')
            warning('Skipping %s: missing field cv_predictions', pred_name);
            continue;
        end
        if ~isfield(predictions, 'all_fullpaths_sorted')
            warning('Skipping %s: missing field all_fullpaths_sorted', pred_name);
            continue;
        end

        all_means       = predictions.cv_predictions(:)';
        image_filenames = predictions.all_fullpaths_sorted;

        n_stim = numel(all_means);
        fprintf('Loaded %d predictions from %s\n', n_stim, pred_name);

        if numel(image_filenames) ~= n_stim
            warning('Skipping %s: number of filenames (%d) does not match number of predictions (%d).', ...
                pred_name, numel(image_filenames), n_stim);
            continue;
        end

        %% ===== EXTRACT STIMULUS NUMBERS FROM FILENAMES =====
        stim_numbers = nan(1, n_stim);
        for i = 1:n_stim
            [~, nameOnly, ~] = fileparts(image_filenames{i});
            num_str = regexp(nameOnly, '\d+', 'match', 'once');
            if ~isempty(num_str)
                stim_numbers(i) = str2double(num_str);
            end
        end

        if any(isnan(stim_numbers))
            warning('%d filenames had no extractable number in %s.', sum(isnan(stim_numbers)), pred_name);
        end

        %% ===== MAP STIMULUS NUMBERS → POSITIONS IN all_means =====
        category_indices = cell(1, n_cat);
        for c = 1:n_cat
            category_indices{c} = find(ismember(stim_numbers, category_stim{c}));
        end

        total_assigned = sum(cellfun(@numel, category_indices));
        fprintf('Assigned %d / %d stimuli to categories\n', total_assigned, n_stim);

        fprintf('\nCategory counts:\n');
        for c = 1:n_cat
            fprintf('  %-10s: %d stimuli (expected %d)\n', ...
                categories{c}, numel(category_indices{c}), numel(category_stim{c}));
        end

        %% ============================================================
        %% 1) RAW ACTIVATIONS
        %% ============================================================
        raw_vals = cellfun(@(idx) all_means(idx(:)), category_indices, 'UniformOutput', false);

        raw_stats = run_metric_pipeline( ...
            raw_vals, categories, bar_colors, clean_name_sp, outDir, ...
            'raw_activation', 'Raw mean activation');

        %% ============================================================
        %% 2) d' VALUES (PER-STIMULUS, POOLED VARIANCE)
        %% ============================================================
        dprime_vals = cell(1, n_cat);
        dprime_rows = {};

        for i = 1:n_cat
            idx_cat = category_indices{i};
            idx_non = setdiff(1:n_stim, idx_cat);

            cat_vals = all_means(idx_cat);
            non_vals = all_means(idx_non);

            mu_non = mean(non_vals, 'omitnan');
            var_cat = var(cat_vals, 'omitnan');
            var_non = var(non_vals, 'omitnan');

            sigma_pooled = sqrt(0.5 * (var_cat + var_non));

            if sigma_pooled == 0 || isnan(sigma_pooled)
                warning('Category %s in %s has pooled SD = 0 or NaN. Setting per-stim d'' to NaN.', ...
                    categories{i}, pred_name);
                vals = nan(size(cat_vals));
            else
                vals = (cat_vals - mu_non) ./ sigma_pooled;
            end

            dprime_vals{i} = vals(:);

            stim_ids_this_cat = stim_numbers(idx_cat);
            pred_vals_this_cat = cat_vals(:);

            for j = 1:numel(vals)
                dprime_rows(end+1, :) = { ...
                    clean_name, ...
                    categories{i}, ...
                    stim_ids_this_cat(j), ...
                    pred_vals_this_cat(j), ...
                    mu_non, ...
                    sigma_pooled, ...
                    vals(j)};
            end
        end

        dprime_stats = run_metric_pipeline( ...
            dprime_vals, categories, bar_colors, clean_name_sp, outDir, ...
            'dprime_selectivity', 'd'' selectivity');

        %% ===== SAVE SINGLE-STIMULUS d' CSV =====
        T_dprime_single = cell2table(dprime_rows, ...
            'VariableNames', {'prediction_file','category','stimulus_number','raw_prediction','mu_non','sigma_pooled','dprime_single_stimulus'});

        writetable(T_dprime_single, fullfile(outDir, 'dprime_selectivity_single_stimulus.csv'));
        fprintf('Saved single-stimulus d'' values → %s\n', ...
            fullfile(outDir, 'dprime_selectivity_single_stimulus.csv'));

        fprintf('\nFinished: %s\n', pred_name);

    catch ME
        warning('Error while processing %s:\n%s', pred_name, ME.message);
    end
end

fprintf('\nAll files processed.\n');

%% =================================================================
%                      HELPER FUNCTIONS
%% =================================================================

function stats_out = run_metric_pipeline(metric_vals, categories, bar_colors, clean_name, outDir, file_tag, ylab)

    [metric_mean, metric_sem, n_vals] = compute_category_summary(metric_vals);

    fprintf('\n=== %s ===\n', strrep(file_tag, '_', ' '));
    for i = 1:numel(categories)
        fprintf('%-10s: mean = %+.3f ± %.3f (SEM), n = %d\n', ...
            categories{i}, metric_mean(i), metric_sem(i), n_vals(i));
    end

    test_vs0 = test_against_zero(metric_vals, numel(categories));

    fprintf('\n=== One-sample t-tests against 0 (Bonferroni corrected) ===\n');
    for i = 1:numel(categories)
        fprintf('%-10s: t(%d) = %.3f, p = %.4f, p_bonf = %.4f\n', ...
            categories{i}, test_vs0.df(i), test_vs0.t(i), test_vs0.p(i), test_vs0.p_bonf(i));
    end

    % Compare highest vs second-highest category distributions
    pair_test = test_highest_vs_second(metric_vals, metric_mean, categories);

    fprintf('\n=== Independent t-test: Highest vs second-highest ===\n');
    fprintf('Highest:        %s (mean = %.3f)\n', categories{pair_test.idx1}, metric_mean(pair_test.idx1));
    fprintf('Second highest: %s (mean = %.3f)\n', categories{pair_test.idx2}, metric_mean(pair_test.idx2));
    fprintf('t(%g) = %.3f, p = %.4f\n', pair_test.df, pair_test.t, pair_test.p);

    save_metric_csvs(outDir, file_tag, categories, n_vals, metric_mean, metric_sem, test_vs0, pair_test);

    plot_metric_barplot(metric_mean, metric_sem, test_vs0.p_bonf, pair_test, ...
        categories, bar_colors, clean_name, ylab, ...
        fullfile(outDir, [file_tag '_barplot.png']));

    stats_out.mean = metric_mean;
    stats_out.sem = metric_sem;
    stats_out.n = n_vals;
    stats_out.test_vs0 = test_vs0;
    stats_out.pair_test = pair_test;
end

function [metric_mean, metric_sem, n_vals] = compute_category_summary(metric_vals)
    n_cat = numel(metric_vals);
    metric_mean = zeros(1, n_cat);
    metric_sem  = zeros(1, n_cat);
    n_vals      = zeros(1, n_cat);

    for i = 1:n_cat
        vals = metric_vals{i};
        vals = vals(~isnan(vals));

        n_vals(i) = numel(vals);

        if isempty(vals)
            metric_mean(i) = NaN;
            metric_sem(i)  = NaN;
        else
            metric_mean(i) = mean(vals, 'omitnan');
            metric_sem(i)  = std(vals, 'omitnan') / sqrt(numel(vals));
        end
    end
end

function test_vs0 = test_against_zero(metric_vals, n_cat)
    p_vals  = nan(1, n_cat);
    t_vals  = nan(1, n_cat);
    df_vals = nan(1, n_cat);

    for i = 1:n_cat
        vals = metric_vals{i};
        vals = vals(~isnan(vals));

        if numel(vals) < 2
            warning('Category %d has fewer than 2 valid values. t-test vs 0 skipped.', i);
            continue;
        end

        [~, p_vals(i), ~, stats] = ttest(vals, 0);
        t_vals(i)  = stats.tstat;
        df_vals(i) = stats.df;
    end

    test_vs0.p = p_vals;
    test_vs0.p_bonf = min(p_vals * n_cat, 1);
    test_vs0.t = t_vals;
    test_vs0.df = df_vals;
end

function pair_test = test_highest_vs_second(metric_vals, metric_mean, categories)
    [~, ord] = sort(metric_mean, 'descend', 'MissingPlacement', 'last');
    idx1 = ord(1);
    idx2 = ord(2);

    vals1 = metric_vals{idx1};
    vals2 = metric_vals{idx2};

    vals1 = vals1(~isnan(vals1));
    vals2 = vals2(~isnan(vals2));

    if numel(vals1) < 2 || numel(vals2) < 2
        warning('Not enough values for highest-vs-second test.');
        p_pair = NaN;
        tstat  = NaN;
        df     = NaN;
        n_used = 0;
    else
        [~, p_pair, ~, stats_pair] = ttest2(vals1, vals2);
        tstat = stats_pair.tstat;
        df    = stats_pair.df;
        n_used = [numel(vals1), numel(vals2)];
    end

    pair_test.idx1 = idx1;
    pair_test.idx2 = idx2;
    pair_test.cat1 = categories{idx1};
    pair_test.cat2 = categories{idx2};
    pair_test.n = n_used;
    pair_test.t = tstat;
    pair_test.df = df;
    pair_test.p = p_pair;
end

function save_metric_csvs(outDir, file_tag, categories, n_vals, metric_mean, metric_sem, test_vs0, pair_test)

    T_metric = table( ...
        categories(:), ...
        n_vals(:), ...
        metric_mean(:), ...
        metric_sem(:), ...
        test_vs0.t(:), ...
        test_vs0.df(:), ...
        test_vs0.p(:), ...
        test_vs0.p_bonf(:), ...
        'VariableNames', {'category','n','mean_value','SEM_value','t_vs_0','df_vs_0','p_vs_0','p_bonf_vs_0'} ...
    );
    writetable(T_metric, fullfile(outDir, [file_tag '_stats.csv']));

    T_pair = table( ...
        string(pair_test.cat1), ...
        string(pair_test.cat2), ...
        metric_mean(pair_test.idx1), ...
        metric_mean(pair_test.idx2), ...
        string(mat2str(pair_test.n)), ...
        pair_test.t, ...
        pair_test.df, ...
        pair_test.p, ...
        'VariableNames', {'highest_category','second_highest_category','mean_highest','mean_second','n_used','t_high_vs_second','df_high_vs_second','p_high_vs_second'} ...
    );
    writetable(T_pair, fullfile(outDir, [file_tag '_pairwise_highest_vs_second.csv']));

    fprintf('Saved stats → %s\n', fullfile(outDir, [file_tag '_stats.csv']));
    fprintf('Saved pairwise stats → %s\n', fullfile(outDir, [file_tag '_pairwise_highest_vs_second.csv']));
end

function plot_metric_barplot(metric_mean, metric_sem, p_bonf, pair_test, categories, bar_colors, clean_name, ylab, outFile)

    n_cat = numel(categories);
    fig = figure('Position', [100 100 700 520]); hold on;

    b = bar(1:n_cat, metric_mean, 'FaceColor', 'flat', 'EdgeColor', 'none');
    for c = 1:n_cat
        b.CData(c,:) = bar_colors(c,:);
    end

    errorbar(1:n_cat, metric_mean, metric_sem, ...
        'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);

    yline(0, 'k--', 'LineWidth', 0.5);

    % fixed limits, intentionally used
    ylim([-1.5 4]);

    yl = ylim;
    y_span = yl(2) - yl(1);
    if y_span == 0
        y_span = 1;
    end
    star_offset = 0.03 * y_span;

    % Per-category stars vs baseline
    for i = 1:n_cat
        if isnan(p_bonf(i)) || isnan(metric_mean(i)) || isnan(metric_sem(i))
            continue;
        end

        if p_bonf(i) < 0.001
            star_str = '***';
        elseif p_bonf(i) < 0.01
            star_str = '**';
        elseif p_bonf(i) < 0.05
            star_str = '*';
        else
            star_str = '';
        end

        if ~isempty(star_str)
            if metric_mean(i) >= 0
                star_y = metric_mean(i) + metric_sem(i) + star_offset;
            else
                star_y = metric_mean(i) - metric_sem(i) - star_offset;
            end

            text(i, star_y, star_str, ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 14, 'FontWeight', 'bold');
        end
    end

    % Highest vs second-highest bracket
    if ~isnan(pair_test.p) && ~isnan(metric_mean(pair_test.idx1)) && ~isnan(metric_mean(pair_test.idx2))
        if pair_test.p < 0.001
            star_pair = '***';
        elseif pair_test.p < 0.01
            star_pair = '**';
        elseif pair_test.p < 0.05
            star_pair = '*';
        else
            star_pair = 'n.s.';
        end

        x1 = pair_test.idx1;
        x2 = pair_test.idx2;

        y_base = max([metric_mean(x1)+metric_sem(x1), metric_mean(x2)+metric_sem(x2)]);
        pad = 0.05 * y_span;
        h   = 0.03 * y_span;
        text_offset = 0.02 * y_span;
        y = y_base + pad;

        plot([x1 x1 x2 x2], [y y+h y+h y], 'k-', 'LineWidth', 1.5);
        text(mean([x1 x2]), y + h + text_offset, star_pair, ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 13, 'FontWeight', 'bold');
    end

    set(gca, ...
        'XTick', 1:n_cat, ...
        'XTickLabel', categories, ...
        'FontSize', 14, ...
        'TickLength', [0 0]);

    ylabel(ylab, 'FontSize', 14);
    title(sprintf('%s: %s', ylab, clean_name), 'FontSize', 14, 'Interpreter', 'none');

    box off;
    hold off;

    exportgraphics(fig, outFile, 'Resolution', 300);
    close(fig);
    fprintf('Saved bar plot → %s\n', outFile);
end