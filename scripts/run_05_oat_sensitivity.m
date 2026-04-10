setup_paths
clear functions
rehash

% -----------------------------
% Settings
% -----------------------------
maxBaseBeats = 500;
maxDrugBeats = 200;

% Main paper analysis: female endocardial T2DM
cases = {
    'endo','female'
};

% Optional supplemental check:
includeMaleCheck = false;
if includeMaleCheck
    cases = [cases; {'endo','male'}];
end

targets = {'none','gkr','inaca','jup'};
scalesByTarget = {
    1.0           % nominal
    [0.8 1.2]     % GKr
    [0.8 1.2]     % INaCa
    [0.8 1.2]     % Jup
};

% Nominal female endocardial threshold from current paper
nominalC_ref = 4.25;

% Local grid around the nominal threshold
localGrid = nominalC_ref * [0.5 0.75 1.0 1.25 1.5 2.0];   % ~2.13 to 8.50 nM

% Fallback grids
lowGrid  = [0.1 0.3 1.0];
highGrid = [10 30 100 300 1000];

summaryRows = [];

for i = 1:size(cases,1)
    celltype = cases{i,1};
    sex      = cases{i,2};

    for j = 1:numel(targets)
        target = targets{j};
        scales = scalesByTarget{j};

        for k = 1:numel(scales)
            scale = scales(k);

            fprintf('\n=== OAT sensitivity: %s | %s | %s | x%.2f ===\n', ...
                celltype, sex, target, scale);

            % -----------------------------
            % Rebuild perturbed baseline from scratch
            % -----------------------------
            [param, X0] = build_base_param(celltype);
            param = apply_sex_layer(param, sex);
            param = apply_t2dm_layer(param, 't2dm', 'none');
            [param, targetLabel] = apply_oat_perturbation(param, target, scale);

            baseOut = pace_to_convergence(X0, param, maxBaseBeats);

            if ~(baseOut.converged && baseOut.isPhys)
                fprintf('  Baseline unusable after perturbation.\n');

                summaryRows = [summaryRows; { ...
                    celltype, sex, char(targetLabel), scale, ...
                    baseOut.converged, baseOut.isPhys, ...
                    baseOut.bm.APD90, baseOut.bm.RMP, ...
                    false, NaN, NaN, '', NaN, NaN, NaN, NaN, NaN}]; %#ok<AGROW>

                outFile = fullfile('results', sprintf( ...
                    'oat_%s_%s_%s_%03d.mat', celltype, sex, lower(char(targetLabel)), round(100*scale)));
                inst = struct();
                save(outFile, 'param', 'baseOut', 'inst');
                continue;
            end

            % -----------------------------
            % Stage 1: local grid around nominal threshold
            % -----------------------------
            inst = find_instability_threshold(baseOut.X0, param, localGrid, maxDrugBeats);

            testedC_all = inst.testedC(:);

            % -----------------------------
            % Stage 2a: if nothing unstable locally, extend upward
            % -----------------------------
            if ~inst.bracketFound
                fullGrid = unique([localGrid highGrid], 'stable');
                inst = find_instability_threshold(baseOut.X0, param, fullGrid, maxDrugBeats);
                testedC_all = inst.testedC(:);
            end

            % -----------------------------
            % Stage 2b: if first local point is already unstable, extend downward
            % -----------------------------
            if inst.bracketFound && ~isnan(inst.C_instab)
                if inst.C_instab <= min(localGrid) + eps
                    fullGrid = unique([lowGrid localGrid], 'stable');
                    inst = find_instability_threshold(baseOut.X0, param, fullGrid, maxDrugBeats);
                    testedC_all = inst.testedC(:);
                end
            end

            % Fold-shift relative to nominal female endo threshold
            foldShift = NaN;
            if ~isnan(inst.C_instab)
                foldShift = inst.C_instab / nominalC_ref;
            end

            % For paper reporting, keep tested range visible
            if isempty(testedC_all)
                testedMin = NaN;
                testedMax = NaN;
                nTests = 0;
            else
                testedMin = min(testedC_all);
                testedMax = max(testedC_all);
                nTests = numel(testedC_all);
            end

            fprintf('  base: converged=%d | isPhys=%d | APD90=%.2f | RMP=%.2f\n', ...
                baseOut.converged, baseOut.isPhys, baseOut.bm.APD90, baseOut.bm.RMP);

            fprintf('  C_instab = %g nM | trigger = %s | fold vs 4.25 nM = %.3f\n', ...
                inst.C_instab, string_or_na(inst.triggerType), foldShift);

            summaryRows = [summaryRows; { ...
                celltype, sex, char(targetLabel), scale, ...
                baseOut.converged, baseOut.isPhys, ...
                baseOut.bm.APD90, baseOut.bm.RMP, ...
                inst.bracketFound, inst.C_instab, inst.beta_instab, ...
                char_or_blank(inst.triggerType), ...
                foldShift, testedMin, testedMax, nTests, ...
                nominalC_ref}]; %#ok<AGROW>

            outFile = fullfile('results', sprintf( ...
                'oat_%s_%s_%s_%03d.mat', celltype, sex, lower(char(targetLabel)), round(100*scale)));
            save(outFile, 'param', 'baseOut', 'inst');
        end
    end
end

oatTable = cell2table(summaryRows, 'VariableNames', { ...
    'CellType','Sex','Target','Scale', ...
    'BaseConverged','BaseIsPhys','BaseAPD90','BaseRMP', ...
    'BracketFound','C_instab_nM','Beta_instab','TriggerType', ...
    'FoldShift_vs_ref','TestedMin_nM','TestedMax_nM','Ntests','ReferenceThreshold_nM'});

disp(oatTable)
writetable(oatTable, fullfile('results', 'oat_sensitivity_summary.csv'));
save(fullfile('results', 'oat_sensitivity_summary.mat'), 'oatTable');

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