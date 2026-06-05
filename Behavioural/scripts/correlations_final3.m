% ============================================================
% SUBJECT RDMs -> correlate to matching original average RDM -> average rhos
% - canonicalizes every loaded RDM into the SAME format first
% - computes subject-wise Spearman correlations
% - runs subject-wise label-shuffle permutation tests
% - applies max-stat correction across conditions
% - plots mean subject rho with individual subject dots
% ============================================================

clear; clc; close all;

%% ============================================================
%% PATHS
%% ============================================================
truePath     = '/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/behavioural/results/RDMs_original/';
sub_basePath = '/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/behavioural/analysis/results_sorted/';
outputDir    = '/Users/lunameidoering/Uni/MASTER/Erasmus/MASTERTHESIS/behavioural/results/correlation_output/';

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% ============================================================
%% SETTINGS
%% ============================================================
nPerm = 10000;
rng(1);

%% ============================================================
%% FIND SUBJECT CONDITION FOLDERS
%% ============================================================
allEntries = dir(sub_basePath);
isValidDir = [allEntries.isdir] & ~ismember({allEntries.name}, {'.','..'});
condFolders = allEntries(isValidDir);

keepMask = false(numel(condFolders),1);
for i = 1:numel(condFolders)
    keepMask(i) = startsWith(condFolders(i).name, 't_') || startsWith(condFolders(i).name, 'g_');
end
condFolders = condFolders(keepMask);

if isempty(condFolders)
    error('No t_* or g_* subject folders found in %s', sub_basePath);
end

fprintf('Found %d subject folders.\n', numel(condFolders));

%% ============================================================
%% LOAD SUBJECT DATA AND COMPUTE OBSERVED SUBJECT CORRELATIONS
%% ============================================================
fprintf('\n============================================================\n');
fprintf('Running subject-level analysis: individual -> original, then average\n');
fprintf('Permutation test: label shuffle within each subject RDM\n');
fprintf('============================================================\n');

subjectResults = table();

subjSqPerCond_map  = struct();
trueVecPerCond_map = struct();

for i = 1:numel(condFolders)

    condName = condFolders(i).name;
    folderPath = fullfile(sub_basePath, condName);

    fprintf('\nProcessing subject folder %s\n', condName);

    parts = split(condName, '_');
    if numel(parts) < 2
        warning('Skipping folder %s: could not parse condition/task.', condName);
        continue;
    end

    stimType = parts{1};
    taskType = parts{2};

    if ~ismember(stimType, {'t','g'})
        warning('Skipping folder %s: stimType is not t or g.', condName);
        continue;
    end

    trueFiles = dir(fullfile(truePath, ['o_' taskType '*.mat']));
    if isempty(trueFiles)
        warning('No matching true file for task %s', taskType);
        continue;
    end

    trueDataLoaded = load(fullfile(truePath, trueFiles(1).name));
    trueRDM_raw = extract_rdm_from_struct(trueDataLoaded);
    [trueSq, trueVec] = canonicalize_rdm(trueRDM_raw);

    subjFiles = dir(fullfile(folderPath, '*.mat'));
    if isempty(subjFiles)
        warning('No .mat files found in %s', folderPath);
        continue;
    end

    fprintf('  Found %d subject files\n', numel(subjFiles));

    sqList = {};

    for s = 1:numel(subjFiles)
        subjFile = subjFiles(s).name;
        subjPath = fullfile(folderPath, subjFile);

        subjLoaded = load(subjPath);
        subjRDM_raw = extract_rdm_from_struct(subjLoaded);
        [subjSq, subjVec] = canonicalize_rdm(subjRDM_raw);

        if numel(subjVec) ~= numel(trueVec)
            warning('Skipping %s/%s: length mismatch (true=%d, subj=%d)', ...
                condName, subjFile, numel(trueVec), numel(subjVec));
            continue;
        end

        if size(subjSq,1) ~= size(trueSq,1)
            warning('Skipping %s/%s: square size mismatch (true=%d, subj=%d)', ...
                condName, subjFile, size(trueSq,1), size(subjSq,1));
            continue;
        end

        rhoSubj = corr(trueVec, subjVec, 'type', 'Spearman', 'rows', 'complete');

        fprintf('    %s: rho = %.4f\n', subjFile, rhoSubj);

        newRow = table({condName}, {make_pretty_label(stimType, taskType)}, ...
            {stimType}, {taskType}, {subjFile}, rhoSubj, ...
            'VariableNames', {'Condition', 'Label', 'StimType', 'TaskType', 'SubjectFile', 'SpearmanRho'});
        subjectResults = [subjectResults; newRow];

        sqList{end+1,1} = subjSq; %#ok<AGROW>
    end

    if ~isempty(sqList)
        subjSqPerCond_map.(condName)  = sqList;
        trueVecPerCond_map.(condName) = trueVec;
    end
end

if isempty(subjectResults)
    error('No valid subject-level correlations were computed.');
end

%% ============================================================
%% AGGREGATE SUBJECT-LEVEL RHOS BY CONDITION
%% ============================================================
uConds = unique(subjectResults.Condition, 'stable');

subjCondNames  = cell(numel(uConds),1);
subjNiceLabels = cell(numel(uConds),1);
subjMeanRhos   = nan(numel(uConds),1);
subjSemRhos    = nan(numel(uConds),1);
subjAllRhos    = cell(numel(uConds),1);
subjN          = nan(numel(uConds),1);
subjPUnc       = nan(numel(uConds),1);

for c = 1:numel(uConds)
    condName = uConds{c};
    mask = strcmp(subjectResults.Condition, condName);
    vals = subjectResults.SpearmanRho(mask);

    subjCondNames{c}  = condName;
    subjNiceLabels{c} = subjectResults.Label{find(mask,1,'first')};
    subjMeanRhos(c)   = mean(vals, 'omitnan');
    subjSemRhos(c)    = std(vals, 'omitnan') / sqrt(sum(~isnan(vals)));
    subjAllRhos{c}    = vals;
    subjN(c)          = sum(~isnan(vals));
end

%% ============================================================
%% SUBJECT-WISE LABEL-SHUFFLE PERMUTATION TEST
%% ============================================================
fprintf('\nRunning subject-wise label-shuffle permutation test...\n');

permMeanMat = nan(nPerm, numel(uConds));

for c = 1:numel(uConds)
    condName = uConds{c};

    if ~isfield(subjSqPerCond_map, condName) || ~isfield(trueVecPerCond_map, condName)
        warning('Missing stored permutation data for %s', condName);
        continue;
    end

    trueVec = trueVecPerCond_map.(condName);
    subjSqList = subjSqPerCond_map.(condName);
    nSubj = numel(subjSqList);

    if nSubj == 0
        continue;
    end

    nItems = size(subjSqList{1}, 1);

    for p = 1:nPerm
        permRhos = nan(nSubj,1);

        for s = 1:nSubj
            subjSq = subjSqList{s};
            permIdx = randperm(nItems);
            permSq = subjSq(permIdx, permIdx);
            permVec = squareform(permSq).';
            permRhos(s) = corr(trueVec, permVec(:), 'type', 'Spearman', 'rows', 'complete');
        end

        permMeanMat(p,c) = mean(permRhos, 'omitnan');
    end

    subjPUnc(c) = (sum(permMeanMat(:,c) >= subjMeanRhos(c)) + 1) / (nPerm + 1);

    fprintf('  %s: mean rho = %.4f, p_unc = %.5f\n', ...
        subjCondNames{c}, subjMeanRhos(c), subjPUnc(c));
end

%% ============================================================
%% MAX-STAT CORRECTION ACROSS CONDITIONS
%% ============================================================
fprintf('\nRunning max-stat correction across conditions...\n');

permMax = max(permMeanMat, [], 2, 'omitnan');
subjPMax = nan(numel(uConds),1);

for c = 1:numel(uConds)
    subjPMax(c) = (sum(permMax >= subjMeanRhos(c)) + 1) / (nPerm + 1);
end

maxStatThresh = prctile(permMax, 95);

fprintf('\n=== Max-stat corrected results ===\n');
for c = 1:numel(uConds)
    fprintf('%s: mean rho = %.4f, p_unc = %.5f, p_max = %.5f\n', ...
        subjCondNames{c}, subjMeanRhos(c), subjPUnc(c), subjPMax(c));
end
fprintf('95th percentile max-null threshold: %.4f\n', maxStatThresh);

%% ============================================================
%% SAVE SUBJECT SUMMARY
%% ============================================================
subjectSummary = table(subjCondNames, subjNiceLabels, subjMeanRhos, subjSemRhos, subjPUnc, subjPMax, subjN, ...
    'VariableNames', {'Condition', 'Label', 'MeanSpearmanRho', 'SEMSpearmanRho', ...
                      'PValueUncorrected', 'PValueMaxStat', 'NSubjects'});

writetable(subjectResults, fullfile(outputDir, 'subject_level_RDM_correlations.csv'));
writetable(subjectSummary, fullfile(outputDir, 'subject_level_RDM_correlations_summary.csv'));

save(fullfile(outputDir, 'subject_level_RDM_correlations_maxstat.mat'), ...
    'subjectResults', 'subjectSummary', 'permMeanMat', 'permMax', 'maxStatThresh');

%% ============================================================
%% SORT BARPLOT FROM HIGHEST TO LOWEST MEAN SUBJECT RHO
%% ============================================================
[plotMeanRhos, sortIdx] = sort(subjMeanRhos, 'descend');
plotSemRhos   = subjSemRhos(sortIdx);
plotLabels    = subjNiceLabels(sortIdx);
plotAllRhos   = subjAllRhos(sortIdx);
plotPMax      = subjPMax(sortIdx);

%% ============================================================
%% BARPLOT: SUBJECT-LEVEL RHOS AVERAGED OVER SUBJECTS
%% ============================================================
figure('Position',[220 220 1000 550]); hold on;

bar(plotMeanRhos, 'FaceColor', [0.65 0.65 0.65], 'EdgeColor', 'none');

errorbar(1:numel(plotMeanRhos), plotMeanRhos, plotSemRhos, ...
    'k.', 'LineWidth', 1.3, 'CapSize', 8);

yline(0, 'k-', 'LineWidth', 0.75);
yline(maxStatThresh, 'r--', '95% max-null threshold', 'LineWidth', 1.5);

DOT_COLOR = [0.1 0.1 0.1];
DOT_SIZE  = 38;

for k = 1:numel(plotAllRhos)
    vals = plotAllRhos{k};
    jitter = (rand(size(vals)) - 0.5) * 0.22;
    scatter(k + jitter, vals, DOT_SIZE, DOT_COLOR, ...
        'filled', 'MarkerFaceAlpha', 0.35);
end

yl = ylim;
ySpan = max(1e-6, yl(2) - yl(1));
starOffset = 0.03 * ySpan;

for k = 1:numel(plotLabels)
    starStr = p_to_stars(plotPMax(k));
    if ~isempty(starStr)
        text(k, plotMeanRhos(k) + plotSemRhos(k) + starOffset, starStr, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 16, ...
            'FontWeight', 'bold');
    end
end

xticks(1:numel(plotLabels));
xticklabels(plotLabels);
xtickangle(35);
ylabel('Mean Subject Spearman correlation to original (\rho)');
box on;
grid on;

saveas(gcf, fullfile(outputDir, 'subject_average_RDM_correlations_barplot_maxstat.png'));

disp(subjectSummary(sortIdx,:));

fprintf('\nSaved subject-level results, summary, max-stat results, and barplot.\n');

%% ============================================================
%% LOCAL FUNCTIONS
%% ============================================================
function RDM = extract_rdm_from_struct(S)
    candidateFields = {'dissimilarityMat', 'estimate_dissimMat_ltv', 'RDM', 'rdm'};
    for f = 1:numel(candidateFields)
        if isfield(S, candidateFields{f})
            RDM = S.(candidateFields{f});
            return;
        end
    end

    fields = fieldnames(S);
    for f = 1:numel(fields)
        x = S.(fields{f});
        if isnumeric(x) && ismatrix(x)
            RDM = x;
            return;
        end
    end

    error('No numeric RDM variable found in loaded .mat file.');
end

function [RDMsq, vec] = canonicalize_rdm(RDM)
    RDMsq = vec_to_square_if_needed(RDM);
    vec = squareform(RDMsq).';
    vec = vec(:);
end

function RDMsq = vec_to_square_if_needed(RDM)
    if ~isvector(RDM)
        if size(RDM,1) ~= size(RDM,2)
            error('RDM matrix must be square.');
        end
        RDMsq = double(RDM);
        RDMsq = (RDMsq + RDMsq.') / 2;
        RDMsq(1:size(RDMsq,1)+1:end) = 0;
        return;
    end

    vec = double(RDM(:));
    m = numel(vec);
    n = (1 + sqrt(1 + 8*m)) / 2;

    if abs(n - round(n)) > 1e-10
        error('Vector length %d does not correspond to a valid upper triangle.', m);
    end
    n = round(n);

    RDMsq = squareform(vec);
    RDMsq = (RDMsq + RDMsq.') / 2;
    RDMsq(1:size(RDMsq,1)+1:end) = 0;
end

function starStr = p_to_stars(p)
    if p < 0.001
        starStr = '***';
    elseif p < 0.01
        starStr = '**';
    elseif p < 0.05
        starStr = '*';
    else
        starStr = '';
    end
end

function label = make_pretty_label(stimType, taskType)
    switch lower(stimType)
        case 'g'
            stimLabel = 'Gatys';
        case 't'
            stimLabel = 'Texform';
        otherwise
            stimLabel = stimType;
    end

    switch lower(taskType)
        case 'vis'
            taskLabel = 'Visual';
        case 'anim'
            taskLabel = 'Animacy';
        case 'siz'
            taskLabel = 'Size';
        otherwise
            taskLabel = taskType;
    end

    label = [stimLabel ' ' taskLabel];
end