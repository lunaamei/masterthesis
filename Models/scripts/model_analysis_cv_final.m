% =========================================================================
% UPDATED FOR CV!
%  - Top/Bottom 50 images
%  - Selectivity Index (SI)
%  - d-prime
%  - Rank-based selectivity
%  - ROC / AUC
%  - Single-stimulus RDMs (predictions + reweighted features)
%  - Category RDMs (predictions + reweighted features)
% =========================================================================

%% ===== CONFIGURATION =====
clear; clc; close all;

addpath('/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/EncodingModels/final/baseline_imagenet_cv/predictions_cv/')
baseOutDir = '/Users/lunameidoering/Desktop/CODE/Models/baseline_imagenet_cv/results/';

pred_name    = "encoding_model_baseline_natural_LOTC_body_left.mat";

image_folder = '/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/image_synthesis/stimuli/natural/';
%image_folder = '/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/image_synthesis/stimuli/gatys_fully_processed/configuration4/';
%image_folder = '/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/image_synthesis/stimuli/texforms/';

outDir = fullfile(baseOutDir, pred_name);
if ~exist(outDir,'dir'); mkdir(outDir); end

N_top = 25;   % stimuli used for rank-based selectivity


%% ===== LOAD MODEL OUTPUT =====
predictions       = load(pred_name);
%all_means         = predictions.final_predictions(:)';
all_means         = predictions.cv_predictions(:)';
image_filenames   = predictions.all_fullpaths_sorted;
reweighted_features = predictions.reweighted_features;   

n_stim = numel(all_means);
fprintf('Loaded %d predictions\n', n_stim);



%% ===== EXTRACT STIMULUS NUMBERS =====
stim_numbers = zeros(1, n_stim);
for i = 1:n_stim
    [~, name, ~] = fileparts(image_filenames{i});
    num_str = regexp(name, '\d+', 'match', 'once');
    stim_numbers(i) = str2double(num_str);
end

if any(isnan(stim_numbers))
    warning('%d filenames had no extractable number!', sum(isnan(stim_numbers)));
end


%% ===== CATEGORY DEFINITIONS =====
bodies_stim = [1:12,   101:113];
hands_stim  = [13:25,  114:125];
tools_stim  = [26:50,  126:150];
mani_stim   = [51:75,  151:175];
nman_stim   = [76:100, 176:200];

categories    = {'bodies','hands','tools','mani','nman'};
category_stim = {bodies_stim, hands_stim, tools_stim, mani_stim, nman_stim};

bar_colors = [
    1.00, 0.50, 0.00;  
    0.30, 0.60, 0.30;  
    0.00, 0.00, 0.50;   
    0.30, 0.75, 0.93;   
    1.00, 0.00, 1.00;  
];


%% ===== MAP STIMULI TO CATEGORIES (unsorted, for d'/ROC) =====
category_indices = cell(1, numel(categories));
for c = 1:numel(categories)
    category_indices{c} = find(ismember(stim_numbers, category_stim{c}));
end

total_assigned = sum(cellfun(@numel, category_indices));
fprintf('Assigned %d / %d stimuli to categories\n', total_assigned, n_stim);

all_assigned       = horzcat(category_indices{:});
unassigned         = setdiff(1:n_stim, all_assigned);
if ~isempty(unassigned)
    warning('%d stimuli could not be assigned to a category.', numel(unassigned));
else
    fprintf('All stimuli successfully assigned.\n');
end


%% ===== SORT BY CATEGORY THEN STIMULUS NUMBER =====
stim_category = strings(n_stim, 1);
cat_order     = zeros(n_stim, 1);

for i = 1:n_stim
    if     ismember(stim_numbers(i), bodies_stim); stim_category(i) = "bodies"; cat_order(i) = 1;
    elseif ismember(stim_numbers(i), hands_stim);  stim_category(i) = "hands";  cat_order(i) = 2;
    elseif ismember(stim_numbers(i), tools_stim);  stim_category(i) = "tools";  cat_order(i) = 3;
    elseif ismember(stim_numbers(i), mani_stim);   stim_category(i) = "mani";   cat_order(i) = 4;
    elseif ismember(stim_numbers(i), nman_stim);   stim_category(i) = "nman";   cat_order(i) = 5;
    end
end

[~, sort_idx]          = sortrows([cat_order(:), stim_numbers(:)], [1 2]);
stim_numbers_sorted    = stim_numbers(sort_idx);
predictions_sorted     = all_means(sort_idx);
features_sorted        = reweighted_features(sort_idx, :);
stim_category_sorted   = stim_category(sort_idx);


%% ===== SORT BY PREDICTION SCORE (for image grids / bar plot) =====
[sorted_means, sorted_idx]    = sort(all_means, 'descend');
sorted_filenames              = image_filenames(sorted_idx);
sorted_stim_numbers           = stim_numbers(sorted_idx);

sorted_categories = cell(1, numel(sorted_stim_numbers));
for i = 1:numel(sorted_stim_numbers)
    sorted_categories{i} = get_category_from_index( ...
        sorted_stim_numbers(i), category_stim, categories);
end

% % %% ================================================================
% % %  TOP 50 IMAGES
% % %% ================================================================
% fprintf('\nPlotting top 50 stimuli...\n')
% 
% figure('Position', [100 100 1400 900]);
% for i = 1:50
%     subplot(5, 10, i);
%     imshow(imread(sorted_filenames{i}));
%     catLabel = get_category_from_index( ...
%         sorted_stim_numbers(i), category_stim, categories);
%     title(sprintf('%s (#%d)\n%.2f', catLabel, sorted_stim_numbers(i), sorted_means(i)), ...
%         'FontSize', 7, 'Interpreter', 'none');
%     axis off;
% end
% sgtitle('Top 50 Stimuli');
% exportgraphics(gcf, fullfile(outDir,'top50_stimuli.png'), 'Resolution', 300);
% 
% 
% %% ================================================================
% %  BOTTOM 50 IMAGES
% %% ================================================================
% fprintf('Plotting bottom 50 stimuli...\n')
% 
% figure('Position', [100 100 1400 900]);
% for i = 1:50
%     rank_idx = n_stim - 50 + i;
%     subplot(5, 10, i);
%     imshow(imread(sorted_filenames{rank_idx}));
%     catLabel = get_category_from_index( ...
%         sorted_stim_numbers(rank_idx), category_stim, categories);
%     title(sprintf('%s (#%d)\n%.2f', catLabel, sorted_stim_numbers(rank_idx), sorted_means(rank_idx)), ...
%         'FontSize', 7, 'Interpreter', 'none');
%     axis off;
% end
% sgtitle('Bottom 50 Stimuli');
% exportgraphics(gcf, fullfile(outDir,'bottom50_stimuli.png'), 'Resolution', 300);
% 
% 
% %% ================================================================
% %  TOP 8 IMAGES
% %% ================================================================
% figure;
% for i = 1:8
%     subplot(2, 4, i);
%     imshow(imread(sorted_filenames{i}));
%     catLabel = get_category_from_index( ...
%         sorted_stim_numbers(i), category_stim, categories);
%     title(sprintf('%s (#%d) | %.2f', catLabel, sorted_stim_numbers(i), sorted_means(i)), ...
%         'FontSize', 8, 'Interpreter', 'none');
%     axis off;
% end
% exportgraphics(gcf, fullfile(outDir,'top8.png'), 'Resolution', 300);


%% ================================================================
%  SORTED STIMULI BAR PLOT
%% ================================================================
figure; hold on;

for i = 1:numel(sorted_means)
    col = [0.5 0.5 0.5];
    for c = 1:numel(category_stim)
        if ismember(sorted_stim_numbers(i), category_stim{c})
            col = bar_colors(c,:);
            break;
        end
    end
    bar(i, sorted_means(i), 'FaceColor', col, 'EdgeColor', 'none');
end

xlim([0, numel(sorted_means)+1]);
ylim([-1 5]);
ylabel('response prediction')
set(gca, 'FontSize', 25, 'TickLength', [0 0], 'XTick', []);

h_leg = gobjects(numel(categories), 1);
for c = 1:numel(categories)
    h_leg(c) = bar(NaN, NaN, 'FaceColor', bar_colors(c,:), 'EdgeColor', 'none');
end
%legend(h_leg, categories, 'Location', 'northeast', 'FontSize', 12);
hold off;
exportgraphics(gcf, fullfile(outDir,'sorted_stimuli_barplot.png'), 'Resolution', 300);


%% ================================================================
%  RANK-BASED SELECTIVITY (top-N)
%% ================================================================
fprintf('\n=== Top-%d Rank-Based Selectivity ===\n', N_top)

topN_stim  = sorted_stim_numbers(1:N_top);
counts     = cellfun(@(x) sum(ismember(topN_stim, x)), category_stim);
proportion = counts(:) ./ N_top * 100;

for i = 1:numel(categories)
    fprintf('  %-10s %d/%d (%.1f%%)\n', categories{i}, counts(i), N_top, proportion(i));
end

writetable( ...
    table(categories(:), counts(:), repmat(N_top,numel(categories),1), proportion(:), ...
          'VariableNames', {'category','count_in_topN','N','percent_in_topN'}), ...
    fullfile(outDir,'rank_based_selectivity.csv'));


%% ================================================================
%  PIE CHART: CATEGORY DISTRIBUTION IN TOP-N (clean)
%% ================================================================
fprintf('Plotting pie chart for top-%d category distribution...\n', N_top)

valid_idx        = counts > 0;
valid_counts     = counts(valid_idx);
valid_prop       = proportion(valid_idx);
valid_colors     = bar_colors(valid_idx,:);

if isempty(valid_counts)
    warning('No nonzero category counts found for top-%d pie chart. Skipping plot.', N_top);
else
    figure('Position', [100 100 600 600]);

    % Percent labels inside slices
    pie_labels = arrayfun(@(x) sprintf('%.0f%%', x), valid_prop, 'UniformOutput', false);
    p = pie(valid_counts, pie_labels);

    % Style slices + enlarge text
    slice_idx = 1;
    for c = 1:numel(valid_counts)
        % Slice color
        p(slice_idx).FaceColor = valid_colors(c,:);
        p(slice_idx).EdgeColor = 'white';
        p(slice_idx).LineWidth = 1;

        % Text styling (inside slice)
        p(slice_idx+1).FontSize = 16;   % <-- make bigger
        p(slice_idx+1).FontWeight = 'bold';
        p(slice_idx+1).Color = 'black';

        slice_idx = slice_idx + 2;
    end

    axis off;  % cleaner look
    title(sprintf('Top-25'), 'FontSize', 20);


    exportgraphics(gcf, ...
        fullfile(outDir, sprintf('top%d_category_piechart.png', N_top)), ...
        'Resolution', 300);
end


%% ================================================================
%  ROC / AUC
% %% ================================================================
% fprintf('\n=== ROC / AUC ===\n')
% 
% figure('Position', [100 100 1200 400]);
% auc_vals = zeros(1, numel(categories));
% 
% for i = 1:numel(categories)
%     labels = zeros(n_stim, 1);
%     labels(category_indices{i}) = 1;
% 
%     [X, Y, ~, AUC] = perfcurve(labels, all_means(:), 1);
%     auc_vals(i) = AUC;
% 
%     subplot(1, numel(categories), i);
%     plot(X, Y, 'LineWidth', 2); hold on;
%     plot([0 1], [0 1], 'k--'); hold off;
%     title(categories{i}, 'FontSize', 10);
%     xlabel('FPR'); ylabel('TPR');
%     axis square; grid on;
%     legend(sprintf('AUC = %.3f', AUC), 'Location', 'southeast');
% 
%     fprintf('  %-10s AUC = %.3f\n', categories{i}, AUC);
% end
% 
% sgtitle('ROC: Category vs. All Others');
% exportgraphics(gcf, fullfile(outDir,'roc_plots.png'), 'Resolution', 300);
% 
% writetable( ...
%     table(categories(:), auc_vals(:), 'VariableNames', {'category','AUC'}), ...
%     fullfile(outDir,'auc_results.csv'));



%% ================================================================
%  FUNCTIONS
%% ================================================================

function catName = get_category_from_index(stim_num, category_stim, categories)
%GET_CATEGORY_FROM_INDEX  Return category name for a stimulus number.
    catName = 'unknown';
    for c = 1:numel(categories)
        if ismember(stim_num, category_stim{c})
            catName = categories{c};
            return
        end
    end
end

