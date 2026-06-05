clear, clc, close all
addpath('C:\Users\dado1\Desktop\texture_fmri\results\multivariate\ds_original\')
indir = 'C:\Users\dado1\Desktop\texture_fmri\results\multivariate\ds_original\';
sub = dir([indir '*_LOTC_hand_left_ds.mat']);  % Get available subjects

subjs = {'sub01', 'sub02', 'sub03', 'sub04', 'sub05', 'sub06', 'sub07', 'sub08', 'sub09', 'sub10', 'sub11', 'sub12', 'sub14'};
numSubjs = numel(subjs);

% Define condition names
conditions.name = {'Bodies', 'Hands', 'Tools', 'Manip', 'Nman'};

% Define colors for each category
colors = [
    1 0.5 0;        % Bodies (Orange)
    1 1 0;          % Hands (Yellow)
    0 0 0.5;        % Tools (Dark Blue)
    0.3 0.75 0.93;  % Mani (Blue)
    1 0 1;  % Scene (Light Blue)
];

% Store activations across all subjects
all_activations = [];

for s = 1:numSubjs
    Bodies_index = [1, 6, 11];
Hands_index = [2, 7, 12];
Tools_index = [3, 8, 13];
Mani_index = [4, 9, 14];
NMan_index = [5, 10, 15];
% Check if subject has an ROI file
    match_idx = find(contains({sub.name}, subjs{s}), 1);  
    if isempty(match_idx)
        warning('No ROI data found for subject %s. Skipping.', subjs{s});
        continue;
    end

    % Load the subject's dataset
    currentfilename = fullfile(indir, sub(match_idx).name);
    currentfile = load(currentfilename);
    
%     % Check if the dataset contains the required variable
%     if ~isfield(currentfile, 'ds') || ~isfield(currentfile.ds, 'samples')
%         warning('Dataset for subject %s does not contain the required field. Skipping.', subjs{s});
%         continue;
%     end

    ds_struct = currentfile.ds.samples;
    ds_invert = ds_struct';  % Transpose for easier indexing

    % Adjust indices for subjects 7, 9, 10 (remove last value)
    if ismember(subjs{s}, {'sub02', 'sub04', 'sub12'})
        Bodies_index = Bodies_index(1:end-1);
        Hands_index = Hands_index(1:end-1);
        Tools_index = Tools_index(1:end-1);
        Mani_index = Mani_index(1:end-1);
        NMan_index = NMan_index(1:end-1);
    end

    % Compute mean activations for each condition
    activations = [
        mean(mean(ds_invert(:, Bodies_index))), ...
        mean(mean(ds_invert(:, Hands_index))), ...
        mean(mean(ds_invert(:, Tools_index))), ...
        mean(mean(ds_invert(:, Mani_index))), ...
        mean(mean(ds_invert(:, NMan_index))), ...
    ];

    % Store activations for group analysis
    all_activations = [all_activations; activations];

%     % Plot the figure for the current subject
%     figure;
%     b = bar(1:5, activations, 'FaceColor', 'flat'); % Using 'flat' to set colors
%     for k = 1:5
%         b.CData(k, :) = colors(k, :); % Apply color to each bar
%     end
% 
%     set(gca, 'XTick', 1:length(conditions.name), 'XTickLabel', conditions.name, 'FontSize', 17);
%     ylim([-0.5 4]);
%     hold on;
% 
%     % Error bars for individual subject
%     subject_sem = std(ds_invert, 0, 1) ./ sqrt(size(ds_invert, 1));
%     errorbar(1:5, activations, subject_sem(1:5), 'k', 'LineStyle', 'none');
% 
%     hold off;
%     ylabel("Average betas");
%     title(['left VOTC - ' subjs{s}], 'FontSize', 25);
end

% Compute group average and SEM
group_mean = mean(all_activations, 1);
group_sem = std(all_activations, 0, 1) / sqrt(size(all_activations, 1));

%% Normalized
% Compute group average and SEM
group_mean = mean(all_activations, 1);
group_sem = std(all_activations, 0, 1) / sqrt(size(all_activations, 1));

% Compute mean across all conditions for each subject
subject_mean_activation = mean(all_activations, 2);  % Mean across conditions for each subject
normalized_activations = all_activations - subject_mean_activation; % Centered data

% Compute new mean and SEM after normalization
group_mean_centered = mean(normalized_activations, 1);
group_sem_centered = std(normalized_activations, 0, 1) / sqrt(size(normalized_activations, 1));

%% Plot Normalized Activations (De-meaned Data)
figure;
b2 = bar(1:5, group_mean_centered, 'FaceColor', 'flat', 'EdgeColor','none');
for k = 1:5
    b2.CData(k, :) = colors(k, :);
end

% % Clean aesthetics
% set(gca, ...
%     'XTick', 1:length(conditions.name), ...
%     'XTickLabel', conditions.name, ...
%     'FontSize', 17, ...
%     'TickLength', [0 0]);  % Removes tick marks

set(gca, 'TickLength', [0 0]);  % Removes tick marks
set(gca, 'Xtick', []);
% set(gca, 'YTickLabel', 'FontSize', 20);

ylim([-1.3 1.1]); % Focused y-axis range
box off;  % Removes the black box outline

hold on;

% Gray vertical error bars + significance annotations
for i = 1:5
    x = i;
    y = group_mean_centered(i);
    err = group_sem_centered(i);
    line([x x], [y-err y+err], 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5); % vertical gray line
    
%     % Add significance text
%     y_offset = 0.05;  % Space above bar for text
%     if group_mean_centered(i) >= 0
%         text(x, y + err + y_offset, sig_labels(i), 'HorizontalAlignment', 'center', 'FontSize', 14);
%     else
%         text(x, y - err - y_offset - 0.05, sig_labels(i), 'HorizontalAlignment', 'center', 'FontSize', 14);
%     end
end
hold off;
ylabel('');
title('');



