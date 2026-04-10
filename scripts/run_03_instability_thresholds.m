setup_paths
clear functions
rehash

load(fullfile('results','baseline_summary.mat'), 'summaryTable');

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

    fprintf('Instability search: %s | %s | %s\n', celltype, sex, disease);

    inFile = fullfile('results', sprintf('baseline_%s_%s_%s.mat', celltype, sex, disease));
    S = load(inFile, 'result');

    X0_ss = S.result.X0_ss;
    baseParam = S.result.param;

    inst = find_instability_threshold(X0_ss, baseParam, concGrid, maxDrugBeats);

    outResult = struct();
    outResult.celltype = celltype;
    outResult.sex = sex;
    outResult.disease = disease;
    outResult.baseParam = baseParam;
    outResult.X0_ss = X0_ss;
    outResult.instability = inst;

    outFile = fullfile('results', sprintf('instability_%s_%s_%s.mat', celltype, sex, disease));
    save(outFile, 'outResult');

    summaryRows = [summaryRows; { ...
        celltype, sex, disease, ...
        inst.bracketFound, inst.C_instab, inst.beta_instab, ...
        inst.lowerC, inst.upperC, char(inst.triggerType), ...
        numel(inst.testedC)}]; %#ok<AGROW>
end

instabilityTable = cell2table(summaryRows, 'VariableNames', { ...
    'CellType','Sex','Disease','BracketFound','C_instab_nM','Beta_instab', ...
    'LowerC_nM','UpperC_nM','TriggerType','Ntests'});

disp(instabilityTable)
writetable(instabilityTable, fullfile('results', 'instability_summary.csv'));
save(fullfile('results', 'instability_summary.mat'), 'instabilityTable');