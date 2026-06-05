clear all
addpath(genpath('/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/MATLAB/CoSMoMVPA'))
config = cosmo_config();

%% Subject IDs
subjs = {'sub11'};

%% ROI suffixes
roi_suffixes = {'_LOTC_body_left', '_LOTC_hand_left', '_LOTC_tool_left'};

%% Paths
base_study_path = '/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/fMRI/multivariate/fMRI_new/fMRI_preprocessed/sub10/texforms/';
base_mask_path = '/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/fMRI/multivariate/fMRI_new/rois/';
base_results_path = '/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/fMRI/multivariate/fMRI_new/ds_texform';


%% Loop over subjects
for s = 1:numel(subjs)
    subj = subjs{s};
    study_path = fullfile(base_study_path, filesep);
    data_fn = fullfile(study_path, 'SPM.mat');

%     if ~exist(data_fn, 'file')
%         warning('SPM.mat not found for subject %s. Skipping.', subj);
%         continue;
%     end

%     % Determine run/trial count dynamically
%     if ismember(subj, {'sub07', 'sub09', 'sub10', 'sub20'})
%         num_runs = 5;
%         num_trials = 40;
%     else
        num_runs = 3;
        num_trials = 15;
%     end

    targets = repmat(1:5, 1, num_runs);
    chunks = floor(((1:num_trials) - 1) / 5) + 1;

    %% Loop over ROIs
    for r = 1:numel(roi_suffixes)
        roi_name = [subj, roi_suffixes{r}];
        mask_fn = fullfile(base_mask_path, [roi_name, '.nii']);

        if ~exist(mask_fn, 'file')
            warning('Mask %s not found. Skipping.', mask_fn);
            continue;
        end

        % Generate dataset
        ds = cosmo_fmri_dataset(data_fn, 'mask', mask_fn, ...
            'targets', targets, 'chunks', chunks);

        cosmo_check_dataset(ds);
        ds = cosmo_remove_useless_data(ds);

        % Save dataset
        out_fn = fullfile(base_results_path, [roi_name, '_ds']);
        save(out_fn, 'ds');

        fprintf('Saved dataset for %s\n', roi_name);
    end
end
