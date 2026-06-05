clear; clc; close all;
%% ================== PARAMETERS ==================

rootPath  = '/Users/lunameidoering/Desktop/CODE/Behavioural/in/results_sorted';
outputDir = '/Users/lunameidoering/Desktop/CODE/Behavioural/out/category_seperability/';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end


task  = 'vis';
conds = {'g','t'};
condLabels = {'Gatys', 'Texform'};

alpha  = 0.05;
colors = [
    0.00 0.45 0.74;   % blue   — Gatys
    0.85 0.33 0.10;   % orange — Texform
];

%% ================== CATEGORY DEFINITIONS ==================
catDefs  = {
    'Body',      1:13;
    'Hand',      14:25;
    'Tool',      26:50;
    'Manip',     51:75;
    'Non-manip', 76:100
};
nCat     = size(catDefs, 1);
catNames = catDefs(:, 1);
nCond    = numel(conds);

%% ================== COMPUTE SEPARABILITY INDEX ==================
sepData = struct();

for c = 1:nCond
    condName = [conds{c} '_' task];
    condPath = fullfile(rootPath, condName);

    if ~exist(condPath, 'dir')
        fprintf('Folder not found: %s\n', condPath);
        continue;
    end

    matFiles = dir(fullfile(condPath, '*.mat'));
    if isempty(matFiles)
        fprintf('No .mat files in: %s\n', condPath);
        continue;
    end

    catVals = cell(nCat, 1);
    for k = 1:nCat, catVals{k} = []; end

    for f = 1:numel(matFiles)
        d = load(fullfile(condPath, matFiles(f).name));

        if isfield(d, 'estimate_dissimMat_ltv')
            RDM = squareform(d.estimate_dissimMat_ltv(:));
        elseif isfield(d, 'dissimilarityMat')
            RDM = d.dissimilarityMat;
        else
            warning('No RDM field in %s — skipping.', matFiles(f).name);
            continue;
        end

        for k = 1:nCat
            idx      = catDefs{k, 2};
            otherIdx = setdiff(1:100, idx);
            W        = RDM(idx, idx);
            within   = mean(W(triu(true(numel(idx)), 1)));
            between  = mean(RDM(idx, otherIdx), 'all');
            catVals{k}(end+1, 1) = between - within;
        end
    end

    sepData.(condName) = catVals;
end

%% ================== ONE-SAMPLE T-TESTS AGAINST ZERO ==================
% Bonferroni across 10 comparisons per condition (5 cats x 2 conds)

nTests = nCat * nCond;   % 10
tStats = struct();

for c = 1:nCond
    condName = [conds{c} '_' task];
    if ~isfield(sepData, condName), continue; end

    for k = 1:nCat
        vals = sepData.(condName){k};
        [~, p, ~, tStat] = ttest(vals, 0);
        tStats.(condName).t(k)     = tStat.tstat;
        tStats.(condName).p_raw(k) = p;
        tStats.(condName).means(k) = mean(vals);
        tStats.(condName).sems(k)  = std(vals) / sqrt(numel(vals));
    end

    p_corr = min(tStats.(condName).p_raw * nTests, 1);
    tStats.(condName).p_corr = p_corr;
    tStats.(condName).sig    = p_corr < alpha;
end

% --- print ---
fprintf('\n=== ONE-SAMPLE T-TESTS: SEPARABILITY vs. 0 (task: %s) ===\n', upper(task));
fprintf('Bonferroni correction across %d tests (alpha = %.3f)\n', nTests, alpha);

for c = 1:nCond
    condName = [conds{c} '_' task];
    if ~isfield(tStats, condName), continue; end

    fprintf('\nCondition: %s\n', condLabels{c});
    fprintf('%-12s | %-8s | %-8s | %-8s | %-10s | %-10s | %-8s\n', ...
        'Category', 'Mean', 'SEM', 't', 'p (raw)', 'p (corr)', 'Result');
    fprintf('%s\n', repmat('-', 1, 76));

    for k = 1:nCat
        result = 'n.s.';
        if tStats.(condName).sig(k), result = 'SIG *'; end
        fprintf('%-12s | %8.3f | %8.3f | %8.3f | %10.4f | %10.4f | %s\n', ...
            catNames{k}, tStats.(condName).means(k), tStats.(condName).sems(k), ...
            tStats.(condName).t(k), tStats.(condName).p_raw(k), ...
            tStats.(condName).p_corr(k), result);
    end
end

%% ================== PAIRED T-TESTS BETWEEN CATEGORIES ==================
% 10 pairwise comparisons per condition, Bonferroni across 10

nPairs      = nCat * (nCat - 1) / 2;
pairedStats = struct();

fprintf('\n=== PAIRED T-TESTS BETWEEN CATEGORIES (task: %s) ===\n', upper(task));
fprintf('Bonferroni correction across %d pairs per condition (alpha = %.3f)\n', nPairs, alpha);

for c = 1:nCond
    condName = [conds{c} '_' task];
    if ~isfield(sepData, condName), continue; end

    pairIdx   = 0;
    p_raw_all = nan(nPairs, 1);
    tVec      = nan(nPairs, 1);
    iVec      = nan(nPairs, 1);
    jVec      = nan(nPairs, 1);

    for i = 1:nCat
        for j = i+1:nCat
            pairIdx = pairIdx + 1;
            [~, p, ~, tStat] = ttest(sepData.(condName){i}, sepData.(condName){j});
            p_raw_all(pairIdx) = p;
            tVec(pairIdx)      = tStat.tstat;
            iVec(pairIdx)      = i;
            jVec(pairIdx)      = j;
        end
    end

    p_corr_all = min(p_raw_all * nPairs, 1);
    sig_all    = p_corr_all < alpha;

    pairedStats.(condName).p_raw  = p_raw_all;
    pairedStats.(condName).p_corr = p_corr_all;
    pairedStats.(condName).sig    = sig_all;
    pairedStats.(condName).t      = tVec;
    pairedStats.(condName).iVec   = iVec;
    pairedStats.(condName).jVec   = jVec;

    fprintf('\nCondition: %s\n', condLabels{c});
    fprintf('%-20s | %-8s | %-10s | %-10s | %-8s\n', ...
        'Comparison', 't', 'p (raw)', 'p (corr)', 'Result');
    fprintf('%s\n', repmat('-', 1, 66));

    for p = 1:nPairs
        compName = sprintf('%s vs %s', catNames{iVec(p)}, catNames{jVec(p)});
        res = 'n.s.';
        if sig_all(p), res = 'SIG *'; end
        fprintf('%-20s | %8.3f | %10.4f | %10.4f | %s\n', ...
            compName, tVec(p), p_raw_all(p), p_corr_all(p), res);
    end
end

%% ================== PLOT: ALL CONDITIONS IN ONE FIGURE ==================
barData = zeros(nCat, nCond);
semData = zeros(nCat, nCond);
sigData = false(nCat, nCond);

for c = 1:nCond
    condName = [conds{c} '_' task];
    if ~isfield(tStats, condName), continue; end
    barData(:, c) = tStats.(condName).means(:);
    semData(:, c) = tStats.(condName).sems(:);
    sigData(:, c) = tStats.(condName).sig(:);
end

figure('Name', ['Separability: ' upper(task)], ...
       'Color', 'w', 'Position', [100 100 700 480]);
hold on;

b = bar(1:nCat, barData, 'grouped', 'EdgeColor', 'none');
for c = 1:nCond
    b(c).FaceColor = colors(c, :);
end

groupWidth = min(0.8, nCond / (nCond + 1.5));
dynOffset  = max(abs(barData(:))) * 0.08;
if dynOffset == 0, dynOffset = 0.02; end

for c = 1:nCond
    xPos = (1:nCat) - groupWidth/2 + (2*c - 1) * groupWidth / (2*nCond);

    errorbar(xPos, barData(:, c), semData(:, c), 'k.', 'LineWidth', 1.2);
    for k = 1:nCat
        if sigData(k, c) && barData(k, c) > 0
            yStar = barData(k, c) + semData(k, c) + dynOffset;
            text(xPos(k), yStar, '*', 'FontSize', 18, ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
    end
end

yline(0, 'k--', 'LineWidth', 1);

yl = ylim;
ylim([yl(1), yl(2) * 1.20]);

ax = gca;
ax.YAxis.Exponent = 0;
set(ax, 'XTick', 1:nCat, 'XTickLabel', catNames, 'FontSize', 12);
ylabel('Separability Index (Between − Within)');
title(['Category Separability: ' upper(task)]);
legend(b, condLabels, 'Location', 'northwest');
grid on;
ax.XGrid = 'off';
box off;

saveas(gcf, fullfile(outputDir, ['separability_' task '_allconds.png']));