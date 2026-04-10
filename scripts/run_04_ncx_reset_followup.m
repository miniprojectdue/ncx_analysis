setup_paths
clear functions
rehash

% Focused NCX-reset follow-up for the main hypothesis:
% compare full T2DM vs T2DM-with-NCX-reset in endocardial male/female cells.

concGrid = [0.1 0.3 1 3 10 30 100 300 1000];
maxBaseBeats = 500;
maxDrugBeats = 200;

cases = {
    'endo','male'
    'endo','female'
};

summaryRows = [];

for i = 1:size(cases,1)
    celltype = cases{i,1};
    sex      = cases{i,2};

    fprintf('\n=== NCX reset follow-up: %s | %s ===\n', celltype, sex);

    % -----------------------------
    % Variant 1: full T2DM
    % -----------------------------
    fullFile = fullfile('results', sprintf('baseline_%s_%s_t2dm.mat', celltype, sex));
    Sfull = load(fullFile, 'result');

    fullParam = Sfull.result.param;
    fullX0ss  = Sfull.result.X0_ss;

    fullInst = find_instability_threshold(fullX0ss, fullParam, concGrid, maxDrugBeats);

    fprintf('Full T2DM:      C_instab = %g nM | trigger = %s\n', ...
        fullInst.C_instab, string_or_na(fullInst.triggerType));

    if isfield(Sfull.result, 'isPhys')
        fullIsPhys = Sfull.result.isPhys;
    else
        fullIsPhys = true;   % fallback for older baseline files
    end
    
    summaryRows = [summaryRows; { ...
        celltype, sex, "FullT2DM", ...
        Sfull.result.converged, fullIsPhys, ...
        Sfull.result.bm.APD90, Sfull.result.bm.RMP, ...
        fullInst.bracketFound, fullInst.C_instab, fullInst.beta_instab, ...
        char_or_blank(fullInst.triggerType)}]; %#ok<AGROW>

    % -----------------------------
    % Variant 2: T2DM with NCX reset to control level
    % -----------------------------
    [paramReset, X0reset] = build_base_param(celltype);
    paramReset = apply_sex_layer(paramReset, sex);
    paramReset = apply_t2dm_layer(paramReset, 't2dm', 'none');

    % Reset only the T2DM NCX multiplier back to control level.
    % Sex-specific baseline NCX embedded in the Holmes male/female model remains intact.
    paramReset.INaCa_Multiplier = 1.0;

    baseReset = pace_to_convergence(X0reset, paramReset, maxBaseBeats);

    fprintf('NCX-reset base: converged=%d | isPhys=%d | APD90=%.2f | RMP=%.2f\n', ...
        baseReset.converged, baseReset.isPhys, baseReset.bm.APD90, baseReset.bm.RMP);

    if baseReset.converged && baseReset.isPhys
        resetInst = find_instability_threshold(baseReset.X0, paramReset, concGrid, maxDrugBeats);

        fprintf('NCX-reset:      C_instab = %g nM | trigger = %s\n', ...
            resetInst.C_instab, string_or_na(resetInst.triggerType));

        summaryRows = [summaryRows; { ...
            celltype, sex, "T2DM_NCXreset", ...
            baseReset.converged, baseReset.isPhys, ...
            baseReset.bm.APD90, baseReset.bm.RMP, ...
            resetInst.bracketFound, resetInst.C_instab, resetInst.beta_instab, ...
            char_or_blank(resetInst.triggerType)}]; %#ok<AGROW>

        % Simple comparison printout
        if ~isnan(fullInst.C_instab) && ~isnan(resetInst.C_instab)
            fprintf('Threshold shift (reset/full) = %.2fx\n', resetInst.C_instab / fullInst.C_instab);
        elseif ~isnan(fullInst.C_instab) && isnan(resetInst.C_instab)
            fprintf('Threshold shift: NCX reset removed instability up to 1000 nM.\n');
        end
    else
        summaryRows = [summaryRows; { ...
            celltype, sex, "T2DM_NCXreset", ...
            baseReset.converged, baseReset.isPhys, ...
            baseReset.bm.APD90, baseReset.bm.RMP, ...
            false, NaN, NaN, ''}]; %#ok<AGROW>

        fprintf('NCX-reset baseline was not usable, threshold search skipped.\n');
    end
end

ncxResetTable = cell2table(summaryRows, 'VariableNames', { ...
    'CellType','Sex','Variant', ...
    'BaseConverged','BaseIsPhys','BaseAPD90','BaseRMP', ...
    'BracketFound','C_instab_nM','Beta_instab','TriggerType'});

disp(ncxResetTable)
writetable(ncxResetTable, fullfile('results', 'ncx_reset_followup.csv'));
save(fullfile('results', 'ncx_reset_followup.mat'), 'ncxResetTable');

% -----------------------------
% local helpers
% -----------------------------
function s = string_or_na(x)
    if isempty(x) || (isstring(x) && strlength(x)==0)
        s = "NA";
    else
        s = string(x);
    end
end

function c = char_or_blank(x)
    if isempty(x) || (isstring(x) && strlength(x)==0)
        c = '';
    else
        c = char(string(x));
    end
end