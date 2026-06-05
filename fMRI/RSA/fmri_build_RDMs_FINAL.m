% ==================================================
% UPDATED FOR NEW SUBJECTS!
% ==================================================

clear; close all; clc;

addpath('/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/fMRI/multivariate/fMRI_new/ds_texform/')
indir = '/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/fMRI/multivariate/fMRI_new/ds_texform/';

outputDir = '/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/fMRI/multivariate/cat_avg_RDMs_13subs/';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% --------------------------------------------------
% Subjects
% --------------------------------------------------
sub = dir([indir '*_LOTC_body_left_ds.mat']);
fig_title = 'texform_LOTC_body.csv';
subjs = {'sub01', 'sub02', 'sub03', 'sub04', 'sub05','sub06', 'sub07', 'sub08', 'sub09', 'sub10', 'sub11','sub12', 'sub14'};
numSubjs = numel(subjs);


% --------------------------------------------------
% Category names
% --------------------------------------------------
catNames = {'Bodies','Hands','Tools','Manip','Non-manip'};
numCats  = numel(catNames);



% --------------------------------------------------
% Collect subject RDMs
% --------------------------------------------------
RDM_all = [];

for s = 1:numSubjs

    fprintf('Processing %s\n', subjs{s});

    % Check if subject has an ROI file
    match_idx = find(contains({sub.name}, subjs{s}), 1);  
    if isempty(match_idx)
        warning('No ROI data found for subject %s. Skipping.', subjs{s});
        continue;
    end

    % % ----------------------------------------------
    % % activate this for original
    % % ----------------------------------------------
    % Bodies_index = [1, 6, 11];
    % Hands_index  = [2, 7, 12];
    % Tools_index  = [3, 8, 13];
    % Mani_index   = [4, 9, 14];
    % NMan_index   = [5, 10, 15];
    % 
    % % subjects with fewer repetitions
    % if ismember(subjs{s}, {'sub02', 'sub04'})
    %     Bodies_index = Bodies_index(1:end-1);
    %     Hands_index  = Hands_index(1:end-1);
    %     Tools_index  = Tools_index(1:end-1);
    %     Mani_index   = Mani_index(1:end-1);
    %     NMan_index   = NMan_index(1:end-1);
    % end

 
    % ----------------------------------------------
    % Load data
    % ----------------------------------------------
    match_idx = find(contains({sub.name}, subjs{s}), 1);
    if isempty(match_idx)
        warning('No file for %s', subjs{s});
        continue;
    end

    data = load(fullfile(indir, sub(match_idx).name));

    % stimulus × voxel matrix
    X = data.ds.samples;   % e.g. 15–20 × nVoxels
    

    
    nTrials = size(X,1);
    
    if nTrials == 20
        Bodies_index = [1, 6, 11, 16];
        Hands_index  = [2, 7, 12, 17];
        Tools_index  = [3, 8, 13, 18];
        Mani_index   = [4, 9, 14, 19];
        NMan_index   = [5, 10, 15, 20];
    
    elseif nTrials == 15
        Bodies_index = [1, 6, 11];
        Hands_index  = [2, 7, 12];
        Tools_index  = [3, 8, 13];
        Mani_index   = [4, 9, 14];
        NMan_index   = [5, 10, 15];
    
    elseif nTrials == 10
        Bodies_index = [1, 6];
        Hands_index  = [2, 7];
        Tools_index  = [3, 8];
        Mani_index   = [4, 9];
        NMan_index   = [5, 10];
    
    else
        warning('Unexpected number of trials for %s: %d. Skipping.', subjs{s}, nTrials);
        continue;
    end

    % Adjust indices for subjects (remove last value)
    if ismember(subjs{s}, {'sub10'})
        Bodies_index = Bodies_index(1:end-1);
        Hands_index = Hands_index(1:end-1);
        Tools_index = Tools_index(1:end-1);
        Mani_index = Mani_index(1:end-1);
        NMan_index = NMan_index(1:end-1);
    end
    
    catIdx = {
        Bodies_index
        Hands_index
        Tools_index
        Mani_index
        NMan_index
    };
    
    X_cat = zeros(numCats, size(X,2));
    
    for c = 1:numCats
        X_cat(c,:) = mean(X(catIdx{c}, :), 1, 'omitnan');
    end
 
    % ----------------------------------------------
    % Category-level RDM
    % ----------------------------------------------
    RDM_cat = squareform(pdist(X_cat, 'correlation'));

    % ----------------------------------------------
    % Convert to table and Save as CSV
    % ----------------------------------------------
    % Convert matrix to table with category names as column headers
    T = array2table(RDM_cat, 'VariableNames', catNames, 'RowNames', catNames);
    
    % Define the output filename for this specific subject
    % e.g., 'sub01_category_RDM.csv'
    filename = fullfile(outputDir, [subjs{s},fig_title]);
    
    % Write the table to CSV (IncludeRowNames ensures headers on the left side)
    writetable(T, filename, 'WriteRowNames', true);
    
    fprintf('Saved RDM to: %s\n', filename);


    % ----------------------------------------------
    % Plot RDM
    % ----------------------------------------------
    fig = figure('Visible','off');
    imagesc(RDM_cat);
    axis square;
    colorbar;

    set(gca, ...
        'XTick', 1:numCats, ...
        'XTickLabel', catNames, ...
        'YTick', 1:numCats, ...
        'YTickLabel', catNames, ...
        'FontSize', 12, ...
        'TickLength',[0 0]);

    xtickangle(45);

    % improve contrast
    clim = prctile(RDM_cat(:), [5 95]);
    caxis(clim);

    % % save
    % outName = fullfile(outputDir, ...
    %     sprintf('%s_fmri_RDM.png', subjs{s}));
    % exportgraphics(fig, outName, 'Resolution', 300);
    % close(fig);

    % ----------------------------------------------
    % Store
    % ----------------------------------------------
    RDM_all(:,:,end+1) = RDM_cat;

end

% --------------------------------------------------
% Average RDM across subjects
% --------------------------------------------------
avgRDM = mean(RDM_all, 3, 'omitnan');
% --------------------------------------------------
% Save group-average category RDM as CSV (with labels)
% --------------------------------------------------
validVarNames = matlab.lang.makeValidName(catNames);
Tavg = array2table(avgRDM, 'VariableNames', validVarNames);
Tavg.Properties.RowNames = catNames;

outCSV = fullfile(outputDir, fig_title);
writetable(Tavg, outCSV, 'WriteRowNames', true);

% --------------------------------------------------
% Plot group-average RDM
% --------------------------------------------------
fig = figure;
imagesc(avgRDM);
axis square;
colorbar;


set(gca, ...
    'XTick', 1:numCats, ...
    'XTickLabel', catNames, ...
    'YTick', 1:numCats, ...
    'YTickLabel', catNames, ...
    'FontSize', 12, ...
    'TickLength',[0 0]);

xtickangle(45);

clim = prctile(avgRDM(:), [5 95]);
caxis(clim);

% exportgraphics(fig, fullfile(outputDir, ...
%     'group_avg_texform_body_fmri_RDM.png'), ...
%     'Resolution', 300);
