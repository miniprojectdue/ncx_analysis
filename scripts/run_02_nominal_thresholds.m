setup_paths
clear functions
rehash

load(fullfile('results','baseline_summary.mat'), 'summaryTable');

% Only proceed with converged + physiological baseline states
validRows = summaryTable.Converged == 1 & summaryTable.IsPhys == 1;

concGrid = [0.1 0.3 1 3 10 30 100 300 1000];
maxDrugBeats = 200;

summaryRows = [];

for i = 1:height(summaryTable)
    if ~validRows(i)
        continue;
    end

    celltype = char(summaryTable.CellType(i));
    sex      = char(summaryTable.Sex(i));
    disease  = char(summaryTable.Disease(i));

    fprintf('Threshold search: %s | %s | %s\n', celltype, sex, disease);

    inFile = fullfile('results', sprintf('baseline_%s_%s_%s.mat', celltype, sex, disease));
    S = load(inFile, 'result');

    X0_ss = S.result.X0_ss;
    baseParam = S.result.param;

    thr = find_ead_threshold(X0_ss, baseParam, concGrid, maxDrugBeats);

    outResult = struct();
    outResult.celltype = celltype;
    outResult.sex = sex;
    outResult.disease = disease;
    outResult.baseParam = baseParam;
    outResult.X0_ss = X0_ss;
    outResult.threshold = thr;

    outFile = fullfile('results', sprintf('threshold_%s_%s_%s.mat', celltype, sex, disease));
    save(outFile, 'outResult');

    summaryRows = [summaryRows; { ...
        celltype, sex, disease, ...
        thr.bracketFound, thr.C_EAD, thr.beta_EAD, ...
        thr.lowerC, thr.upperC, ...
        numel(thr.testedC)}]; %#ok<AGROW>
end

thresholdTable = cell2table(summaryRows, 'VariableNames', { ...
    'CellType','Sex','Disease','BracketFound','C_EAD_nM','Beta_EAD', ...
    'LowerC_nM','UpperC_nM','Ntests'});

disp(thresholdTable)
writetable(thresholdTable, fullfile('results', 'threshold_summary.csv'));
save(fullfile('results', 'threshold_summary.mat'), 'thresholdTable');