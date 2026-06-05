%% ============================================================
%  Build & visualize RDMs for all .mat files in subfolders,
%  Also computes and saves an average RDM per subfolder.
% =============================================================
clear; close all; clc;

%% ---- Paths ----
rootPath = '/Users/lunameidoering/Desktop/CODE/Behavioural/in/results_sorted';
outputPath = '/Users/lunameidoering/Desktop/CODE/Behavioural/out/behavioural_RDMs';
if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end

%% ---- Get subfolders ----
subFolders = dir(rootPath);
subFolders = subFolders([subFolders.isdir] & ~ismember({subFolders.name},{'.','..'}));

%% ---- Loop over subfolders ----
for sf = 1:numel(subFolders)
    subName = subFolders(sf).name;
    subPath = fullfile(rootPath, subName);
    fprintf('\nProcessing folder: %s\n', subName);

    matFiles = dir(fullfile(subPath, '*.mat'));
    if isempty(matFiles)
        fprintf('No .mat files found in %s\n', subName);
        continue;
    end

    RDMs = {};          % full symmetric square RDMs
    dissVecs = {};      % vector-form dissimilarities in squareform order

    %% ---- Loop through each .mat file ----
    for f = 1:numel(matFiles)
        filePath = fullfile(subPath, matFiles(f).name);
        fprintf('   Loading %s\n', matFiles(f).name);

        data = load(filePath);

        % --- Find the dissimilarity representation ---
        if isfield(data, 'estimate_dissimMat_ltv')
            rawDiss = data.estimate_dissimMat_ltv;
        elseif isfield(data, 'dissimilarityMat')
            rawDiss = data.dissimilarityMat;
        else
            warning('No dissimilarity data found in %s', matFiles(f).name);
            continue;
        end

        % --- Convert to canonical full symmetric square RDM ---
        [RDM, dissVec] = canonicalize_rdm(rawDiss);

        % --- Store ---
        RDMs{end+1} = RDM; %#ok<AGROW>
        dissVecs{end+1} = dissVec; %#ok<AGROW>

        %% ---- Save each RDM in consistent format ----
        [~, baseName] = fileparts(matFiles(f).name);
        prefix = [subName '_'];

        % Save .mat in vector format consistent with subject files
        dissimilarityMat = dissVec; %#ok<NASGU>
        RDM_square = RDM; %#ok<NASGU>
        saveNameMat = fullfile(outputPath, [prefix baseName '_RDM.mat']);
        save(saveNameMat, 'dissimilarityMat', 'RDM_square');

        % Plot and save PNG of full square RDM
        fig = figure('Visible','off');
        imagesc(RDM);
        axis square; colorbar;
        title(sprintf('RDM: %s%s', prefix, baseName), 'Interpreter', 'none');
        exportgraphics(fig, fullfile(outputPath, [prefix baseName '_RDM.png']), 'Resolution', 300);
        close(fig);
    end

    %% ---- Compute and save average RDM for this subfolder ----
    if ~isempty(RDMs)
        n = size(RDMs{1},1);
        avgRDM = zeros(n);

        for i = 1:numel(RDMs)
            if size(RDMs{i},1) ~= n || size(RDMs{i},2) ~= n
                error('RDM size mismatch in folder %s.', subName);
            end
            avgRDM = avgRDM + RDMs{i};
        end

        avgRDM = avgRDM / numel(RDMs);

        % Enforce symmetry and zero diagonal after averaging
        avgRDM = (avgRDM + avgRDM.') / 2;
        avgRDM(1:n+1:end) = 0;

        % Convert averaged full square RDM to vector in squareform order
        avgVec = squareform(avgRDM).';
        avgVec = avgVec(:);

        prefix = [subName '_'];

        % Save average in SAME format as subject files
        dissimilarityMat = avgVec; %#ok<NASGU>
        RDM_square = avgRDM; %#ok<NASGU>
        save(fullfile(outputPath, [prefix 'average_RDM.mat']), 'dissimilarityMat', 'RDM_square');

        % Save average RDM as PNG
        fig = figure('Visible','off');
        imagesc(avgRDM);
        axis square; colorbar;
        title(sprintf('Average RDM: %s', subName), 'Interpreter', 'none');
        exportgraphics(fig, fullfile(outputPath, [prefix 'average_RDM.png']), 'Resolution', 300);
        close(fig);

        fprintf('  Saved %d RDMs and 1 average RDM for %s\n', numel(RDMs), subName);
    end
end

fprintf('\nAll done! RDMs saved in: %s\n', outputPath);

%% ---- Local function ----
function [RDMsq, dissVec] = canonicalize_rdm(rawDiss)
    % Convert either:
    % 1) vector-form dissimilarity data
    % 2) square RDM
    % into:
    % - full symmetric square RDM
    % - vector in squareform order

    if isvector(rawDiss)
        dissVec = double(rawDiss(:));
        m = numel(dissVec);
        n = (1 + sqrt(1 + 8*m)) / 2;

        if abs(n - round(n)) > 1e-10
            error('Vector length %d does not correspond to a valid upper triangle.', m);
        end
        n = round(n);

        RDMsq = squareform(dissVec);
        RDMsq = (RDMsq + RDMsq.') / 2;
        RDMsq(1:n+1:end) = 0;

        % re-export vector to guarantee standard order
        dissVec = squareform(RDMsq).';
        dissVec = dissVec(:);

    else
        if size(rawDiss,1) ~= size(rawDiss,2)
            error('Matrix RDM must be square.');
        end

        RDMsq = double(rawDiss);
        RDMsq = (RDMsq + RDMsq.') / 2;
        RDMsq(1:size(RDMsq,1)+1:end) = 0;

        dissVec = squareform(RDMsq).';
        dissVec = dissVec(:);
    end
end