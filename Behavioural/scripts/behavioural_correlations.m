% ============================================================
% SUBJECT RDMs -> correlate to matching original average RDM
% - computes subject-wise Pearson correlations
% - one-sample t-tests against 0
% - paired t-tests (Gatys vs Texform)
% - Bonferroni correction
% - Plot 1: sorted bars with one-sample significance
% - Plot 2: grouped bars with paired-test brackets
% ============================================================

clear; clc; close all;

%% ============================================================
%% PATHS
%% ============================================================
truePath     = '/Users/lunameidoering/Desktop/CODE/Behavioural/in/RDMs_original/';
sub_basePath = '/Users/lunameidoering/Desktop/CODE/Behavioural/in/RDMs_gt/';
outputDir    = '/Users/lunameidoering/Desktop/CODE/Behavioural/out/correlations/';

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% ============================================================
%% FIND SUBJECT CONDITION FOLDERS
%% ============================================================
allEntries = dir(sub_basePath);

isValidDir = [allEntries.isdir] & ...
    ~ismember({allEntries.name}, {'.','..'});

condFolders = allEntries(isValidDir);

keepMask = false(numel(condFolders),1);

for i = 1:numel(condFolders)
    keepMask(i) = startsWith(condFolders(i).name, 't_') || ...
                  startsWith(condFolders(i).name, 'g_');
end

condFolders = condFolders(keepMask);

if isempty(condFolders)
    error('No valid condition folders found.');
end

fprintf('Found %d condition folders.\n', numel(condFolders));

%% ============================================================
%% LOAD SUBJECT DATA
%% ============================================================
subjectResults = table();

for i = 1:numel(condFolders)

    condName = condFolders(i).name;
    folderPath = fullfile(sub_basePath, condName);

    fprintf('\nProcessing %s\n', condName);

    parts = split(condName, '_');

    if numel(parts) < 2
        continue;
    end

    stimType = parts{1};
    taskType = parts{2};

    %% load matching original RDM
    trueFiles = dir(fullfile(truePath, ['o_' taskType '*.mat']));

    if isempty(trueFiles)
        warning('No original RDM found for task %s', taskType);
        continue;
    end

    trueLoaded = load(fullfile(truePath, trueFiles(1).name));

    trueRDM_raw = extract_rdm_from_struct(trueLoaded);

    [~, trueVec] = canonicalize_rdm(trueRDM_raw);

    %% subject files
    subjFiles = dir(fullfile(folderPath, '*.mat'));

    if isempty(subjFiles)
        warning('No subject files in %s', folderPath);
        continue;
    end

    for s = 1:numel(subjFiles)

        subjFile = subjFiles(s).name;
        subjPath = fullfile(folderPath, subjFile);

        subjLoaded = load(subjPath);

        subjRDM_raw = extract_rdm_from_struct(subjLoaded);

        [~, subjVec] = canonicalize_rdm(subjRDM_raw);

        if numel(subjVec) ~= numel(trueVec)
            warning('Skipping %s (vector length mismatch)', subjFile);
            continue;
        end

        rhoSubj = corr(trueVec, subjVec, ...
            'type', 'Pearson', ...
            'rows', 'complete');

        fprintf('  %s: rho = %.4f\n', subjFile, rhoSubj);

        newRow = table( ...
            {condName}, ...
            {make_pretty_label(stimType, taskType)}, ...
            {stimType}, ...
            {taskType}, ...
            {subjFile}, ...
            rhoSubj, ...
            'VariableNames', ...
            {'Condition', ...
             'Label', ...
             'StimType', ...
             'TaskType', ...
             'SubjectFile', ...
             'PearsonRho'});

        subjectResults = [subjectResults; newRow];
    end
end

if isempty(subjectResults)
    error('No valid subject correlations computed.');
end

%% ============================================================
%% CONDITION SUMMARY + ONE-SAMPLE T-TESTS
%% ============================================================
uConds = unique(subjectResults.Condition, 'stable');

nConditions = numel(uConds);

subjCondNames  = cell(nConditions,1);
subjNiceLabels = cell(nConditions,1);

subjMeanRhos = nan(nConditions,1);
subjSemRhos  = nan(nConditions,1);

subjAllRhos = cell(nConditions,1);

subjN = nan(nConditions,1);

tStats_1samp = nan(nConditions,1);
df_1samp     = nan(nConditions,1);
pRaw_1samp   = nan(nConditions,1);

subjStimTypes = cell(nConditions,1);
subjTaskTypes = cell(nConditions,1);

fprintf('\n====================================================\n');
fprintf('ONE-SAMPLE T-TESTS AGAINST 0\n');
fprintf('====================================================\n');

for c = 1:nConditions

    condName = uConds{c};

    mask = strcmp(subjectResults.Condition, condName);

    vals = subjectResults.PearsonRho(mask);

    subjCondNames{c}  = condName;
    subjNiceLabels{c} = subjectResults.Label{find(mask,1,'first')};

    subjMeanRhos(c) = mean(vals, 'omitnan');

    subjSemRhos(c) = std(vals, 'omitnan') ./ ...
        sqrt(sum(~isnan(vals)));

    subjAllRhos{c} = vals;

    subjN(c) = sum(~isnan(vals));

    subjStimTypes{c} = ...
        subjectResults.StimType{find(mask,1,'first')};

    subjTaskTypes{c} = ...
        subjectResults.TaskType{find(mask,1,'first')};

    % use all available data; ttest handles NaNs via 'omitnan' equivalent
    vals = vals(~isnan(vals));

    [~, pRaw_1samp(c), ~, stats] = ttest(vals, 0);

    tStats_1samp(c) = stats.tstat;
    df_1samp(c)     = stats.df;

    fprintf('%s: t(%d)=%.3f, p_raw=%.5f\n', ...
        condName, stats.df, stats.tstat, pRaw_1samp(c));
end

%% ============================================================
%% BONFERRONI CORRECTION (ONE-SAMPLE)
%% ============================================================
% Correct for nStim x nTasks = 2 x 3 = 6 comparisons
%nBonfComparisons_OneSample = numel(unique(subjStimTypes)) * numel(unique(subjTaskTypes));
nBonfComparisons_OneSample = 6;

pBonf_1samp = min(1, ...
    pRaw_1samp .* nBonfComparisons_OneSample);

fprintf('\nBonferroni correction across %d comparisons (nStim x nTasks)\n', ...
    nBonfComparisons_OneSample);

for c = 1:nConditions

    fprintf('%s: p_raw=%.5f | p_bonf=%.5f\n', ...
        subjCondNames{c}, ...
        pRaw_1samp(c), ...
        pBonf_1samp(c));
end

%% ============================================================
%% SAVE ONE-SAMPLE RESULTS
%% ============================================================
subjectSummary = table( ...
    subjCondNames, ...
    subjNiceLabels, ...
    subjMeanRhos, ...
    subjSemRhos, ...
    tStats_1samp, ...
    df_1samp, ...
    pRaw_1samp, ...
    pBonf_1samp, ...
    repmat(nBonfComparisons_OneSample, nConditions, 1), ...
    subjN, ...
    'VariableNames', ...
    {'Condition', ...
     'Label', ...
     'MeanPearsonRho', ...
     'SEMPearsonRho', ...
     'tStatistic', ...
     'df', ...
     'PValueRaw', ...
     'PValueBonferroni', ...
     'NBonferroniComparisons', ...
     'NSubjects'});

writetable(subjectResults, ...
    fullfile(outputDir, ...
    'subject_level_RDM_correlations.csv'));

writetable(subjectSummary, ...
    fullfile(outputDir, ...
    'subject_level_RDM_correlations_summary.csv'));

%% ============================================================
%% PLOT 1
%% ============================================================
[plotMeanRhos, sortIdx] = ...
    sort(subjMeanRhos, 'descend');

plotSemRhos = subjSemRhos(sortIdx);
plotLabels  = subjNiceLabels(sortIdx);
plotAllRhos = subjAllRhos(sortIdx);
plotPBonf   = pBonf_1samp(sortIdx);

figure('Position',[220 220 1000 550]); hold on;

bar(plotMeanRhos, ...
    'FaceColor', [0.65 0.65 0.65], ...
    'EdgeColor', 'none');

errorbar(1:numel(plotMeanRhos), ...
    plotMeanRhos, ...
    plotSemRhos, ...
    'k.', ...
    'LineWidth', 1.3, ...
    'CapSize', 8);

yline(0, 'k-', 'LineWidth', 0.75);

for k = 1:numel(plotAllRhos)

    vals = plotAllRhos{k};

    jitter = (rand(size(vals)) - 0.5) * 0.22;

    scatter(k + jitter, vals, 38, ...
        [0.1 0.1 0.1], ...
        'filled', ...
        'MarkerFaceAlpha', 0.35);
end

yl = ylim;
ySpan = yl(2) - yl(1);

for k = 1:numel(plotLabels)

    starStr = p_to_stars(plotPBonf(k));

    if isempty(starStr)
        continue;
    end

    text(k, ...
        plotMeanRhos(k) + ...
        plotSemRhos(k) + ...
        0.03*ySpan, ...
        starStr, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 16, ...
        'FontWeight', 'bold');
end

xticks(1:numel(plotLabels));
xticklabels(plotLabels);

xtickangle(35);

ylabel('Mean Subject Pearson correlation to original (\rho)');

title(sprintf(['One-sample t-tests vs 0 | ', ...
    'Bonferroni corrected (\\alpha = 0.05 / %d)'], ...
    nBonfComparisons_OneSample));

box on;
grid on;

saveas(gcf, ...
    fullfile(outputDir, ...
    'subject_average_RDM_correlations_barplot.png'));

fprintf('\nSaved Plot 1\n');

%% ============================================================
%% PLOT 2: GROUPED BY TASK
%% ============================================================

uTasks = unique(subjTaskTypes, 'stable');

preferredOrder = {'vis','anim','siz'};

orderedTasks = {};

for k = 1:numel(preferredOrder)

    if any(strcmp(uTasks, preferredOrder{k}))
        orderedTasks{end+1} = preferredOrder{k};
    end
end

for k = 1:numel(uTasks)

    if ~any(strcmp(orderedTasks, uTasks{k}))
        orderedTasks{end+1} = uTasks{k};
    end
end

uTasks = orderedTasks(:);

nTasks = numel(uTasks);

COLOR_G = [0.894 0.447 0.102];
COLOR_T = [0.173 0.482 0.714];

BAR_WIDTH = 0.32;

taskMean   = nan(nTasks,2);
taskSEM    = nan(nTasks,2);
taskAllRho = cell(nTasks,2);

tStatsPair = nan(nTasks,1);
pRawPair   = nan(nTasks,1);

fprintf('\n====================================================\n');
fprintf('PAIRED T-TESTS: GATYS vs TEXFORM\n');
fprintf('====================================================\n');

for ti = 1:nTasks

    task = uTasks{ti};

    %% collect data
    for si = 1:2

        if si == 1
            stim = 'g';
        else
            stim = 't';
        end

        cidx = find(strcmp(subjCondNames, ...
            [stim '_' task]), 1);

        if isempty(cidx)
            continue;
        end

        taskMean(ti,si) = subjMeanRhos(cidx);
        taskSEM(ti,si)  = subjSemRhos(cidx);

        taskAllRho{ti,si} = subjAllRhos{cidx};
    end

    %% paired t-test
    maskG = strcmp(subjectResults.StimType,'g') & ...
            strcmp(subjectResults.TaskType,task);

    maskT = strcmp(subjectResults.StimType,'t') & ...
            strcmp(subjectResults.TaskType,task);

    filesG = subjectResults.SubjectFile(maskG);
    filesT = subjectResults.SubjectFile(maskT);

    rhosG = subjectResults.PearsonRho(maskG);
    rhosT = subjectResults.PearsonRho(maskT);

    % Match subjects by ID = everything before the first underscore
    get_subj_id = @(f) strtok(f, '_');
    baseG = cellfun(get_subj_id, filesG, 'UniformOutput', false);
    baseT = cellfun(get_subj_id, filesT, 'UniformOutput', false);

    [~, iaG, iaT] = intersect(baseG, baseT);

    if numel(iaG) < 2
        warning('%s: not enough matched subjects', task);
        continue;
    end

    [~, pRawPair(ti), ~, stats] = ...
        ttest(rhosG(iaG), rhosT(iaT));

    tStatsPair(ti) = stats.tstat;

    fprintf('%s: t(%d)=%.3f, p_raw=%.5f\n', ...
        task, stats.df, stats.tstat, pRawPair(ti));
end

%% ============================================================
%% BONFERRONI CORRECTION (PAIRED)
%% ============================================================
nValidTests = sum(~isnan(pRawPair));

pBonfPair = min(1, pRawPair .* nValidTests);

fprintf('\nBonferroni correction across %d comparisons\n', ...
    nValidTests);

for ti = 1:nTasks

    fprintf('%s: p_raw=%.5f | p_bonf=%.5f\n', ...
        uTasks{ti}, ...
        pRawPair(ti), ...
        pBonfPair(ti));
end

%% ============================================================
%% SAVE PAIRED T-TEST RESULTS
%% ============================================================
ttestResults = table( ...
    uTasks, ...
    tStatsPair, ...
    pRawPair, ...
    pBonfPair, ...
    repmat(nValidTests, nTasks, 1), ...
    'VariableNames', ...
    {'Task', ...
     'tStatistic', ...
     'PValueRaw', ...
     'PValueBonferroni', ...
     'NBonferroniComparisons'});

writetable(ttestResults, ...
    fullfile(outputDir, ...
    'paired_ttest_stim_within_task.csv'));

disp(ttestResults);

%% ============================================================
%% GROUPED BARPLOT
%% ============================================================
figure('Position',[260 260 800 580]); hold on;

xCentres = 1:nTasks;

xG = xCentres - BAR_WIDTH/2;
xT = xCentres + BAR_WIDTH/2;

bG = bar(xG, taskMean(:,1), BAR_WIDTH, ...
    'FaceColor', COLOR_G, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.85);

bT = bar(xT, taskMean(:,2), BAR_WIDTH, ...
    'FaceColor', COLOR_T, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.85);

errorbar(xG, taskMean(:,1), taskSEM(:,1), ...
    'k.', 'LineWidth', 1.2, 'CapSize', 6);

errorbar(xT, taskMean(:,2), taskSEM(:,2), ...
    'k.', 'LineWidth', 1.2, 'CapSize', 6);

yline(0, 'k-', 'LineWidth', 0.75);

%% one-sample significance stars above each bar
yl_pre = ylim;
ySpan_pre = yl_pre(2) - yl_pre(1);
permStarOffset = 0.03 * ySpan_pre;

for ti = 1:nTasks
    task = uTasks{ti};
    for si = 1:2
        if si == 1; stim = 'g'; xc = xG(ti);
        else;       stim = 't'; xc = xT(ti);
        end
        cidx = find(strcmp(subjCondNames, [stim '_' task]), 1);
        if isempty(cidx); continue; end
        starStr = p_to_stars(pBonf_1samp(cidx));
        if isempty(starStr); continue; end
        ystar = taskMean(ti,si) + taskSEM(ti,si) + permStarOffset;
        text(xc, ystar, starStr, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment',   'bottom', ...
            'FontSize', 12, 'FontWeight', 'bold', ...
            'Color', [0.2 0.2 0.2]);
    end
end

%% subject dots
for ti = 1:nTasks

    for si = 1:2

        if si == 1
            xc = xG(ti);
            dc = COLOR_G;
        else
            xc = xT(ti);
            dc = COLOR_T;
        end

        vals = taskAllRho{ti,si};

        if isempty(vals)
            continue;
        end

        jitter = (rand(size(vals)) - 0.5) * 0.15;

        scatter(xc + jitter, vals, 32, dc, ...
            'filled', ...
            'MarkerFaceAlpha', 0.35, ...
            'MarkerEdgeColor', 'none');
    end
end

%% ============================================================
%% SIGNIFICANCE BRACKETS
%% ============================================================
yl = ylim;
ySpan = yl(2) - yl(1);

bracketLift   = 0.06 * ySpan;
bracketHeight = 0.03 * ySpan;
starLift      = 0.015 * ySpan;

maxBracketY = yl(2);

for ti = 1:nTasks

    starStr = p_to_stars(pBonfPair(ti));

    if isempty(starStr)
        continue;
    end

    yBase = max( ...
        taskMean(ti,1) + taskSEM(ti,1), ...
        taskMean(ti,2) + taskSEM(ti,2)) ...
        + bracketLift;

    yTip = yBase + bracketHeight;

    plot([xG(ti) xG(ti)], [yBase yTip], 'k-', 'LineWidth', 1.3);
    plot([xG(ti) xT(ti)], [yTip  yTip], 'k-', 'LineWidth', 1.3);
    plot([xT(ti) xT(ti)], [yTip yBase], 'k-', 'LineWidth', 1.3);

    text(mean([xG(ti) xT(ti)]), ...
        yTip + starLift, ...
        starStr, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 16, ...
        'FontWeight', 'bold');

    maxBracketY = max(maxBracketY, yTip + 3*starLift);
end

ylim([yl(1), maxBracketY + 0.05*ySpan]);

xticks(xCentres);

xticklabels(cellfun(@task_pretty_label, ...
    uTasks, ...
    'UniformOutput', false));

ylabel('Mean Subject Pearson correlation to natural RDM');

title(sprintf(['Gatys vs Texform | ', ...
    'paired t-tests, Bonferroni corrected ', ...
    '(\alpha = 0.05 / %d)'], ...
    nValidTests));

legend([bG bT], ...
    {'Gatys','Texform'}, ...
    'Location','best', ...
    'Box','off');

box on;
grid on;

saveas(gcf, ...
    fullfile(outputDir, ...
    'subject_RDM_correlations_grouped_ttest.png'));

fprintf('\nSaved Plot 2\n');

%% ============================================================
%% LOCAL FUNCTIONS
%% ============================================================

function RDM = extract_rdm_from_struct(S)

candidateFields = { ...
    'dissimilarityMat', ...
    'estimate_dissimMat_ltv', ...
    'RDM', ...
    'rdm'};

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

error('No numeric RDM found.');
end

function [RDMsq, vec] = canonicalize_rdm(RDM)
RDMsq = vec_to_square_if_needed(RDM);
vec = squareform(RDMsq).';
vec = vec(:);
end

function RDMsq = vec_to_square_if_needed(RDM)

if ~isvector(RDM)
    if size(RDM,1) ~= size(RDM,2)
        error('RDM must be square.');
    end
    RDMsq = double(RDM);
    RDMsq = (RDMsq + RDMsq.') ./ 2;
    RDMsq(1:size(RDMsq,1)+1:end) = 0;
    return;
end

vec = double(RDM(:));
RDMsq = squareform(vec);
RDMsq = (RDMsq + RDMsq.') ./ 2;
RDMsq(1:size(RDMsq,1)+1:end) = 0;
end

function starStr = p_to_stars(p)
if isnan(p)
    starStr = '';
elseif p < 0.001
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
    case 'g';  stimLabel = 'Gatys';
    case 't';  stimLabel = 'Texform';
    otherwise; stimLabel = stimType;
end
switch lower(taskType)
    case 'vis';  taskLabel = 'General';
    case 'anim'; taskLabel = 'Animacy';
    case 'siz';  taskLabel = 'Size';
    otherwise;   taskLabel = taskType;
end
label = [stimLabel ' ' taskLabel];
end

function label = task_pretty_label(taskType)
switch lower(taskType)
    case 'vis';  label = 'General';
    case 'anim'; label = 'Animacy';
    case 'siz';  label = 'Size';
    otherwise;   label = taskType;
end
end