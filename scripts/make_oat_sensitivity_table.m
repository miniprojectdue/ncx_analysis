setup_paths
clear functions
rehash

Tin = readtable(fullfile('results','oat_sensitivity_summary.csv'), 'TextType', 'string');

% Keep main paper case only
idx = strcmpi(Tin.CellType, 'endo') & strcmpi(Tin.Sex, 'female');
T = Tin(idx,:);

% Stable ordering
orderLabels = [
    "Nominal_x1.00"
    "GKr_x0.80"
    "GKr_x1.20"
    "INaCa_x0.80"
    "INaCa_x1.20"
    "Jup_x0.80"
    "Jup_x1.20"
];

rowKey = strings(height(T),1);
for i = 1:height(T)
    rowKey(i) = sprintf('%s_x%.2f', string(T.Target(i)), double(T.Scale(i)));
end

[~, loc] = ismember(orderLabels, rowKey);
loc = loc(loc > 0);
T = T(loc,:);

Perturbation = strings(height(T),1);
C_instab = strings(height(T),1);
FoldShift = strings(height(T),1);
Trigger = strings(height(T),1);
BaseStatus = strings(height(T),1);

for i = 1:height(T)
    if strcmpi(T.Target(i), "Nominal")
        Perturbation(i) = "Nominal";
    else
        signStr = "+20%";
        if double(T.Scale(i)) < 1
            signStr = "-20%";
        end
        Perturbation(i) = sprintf('%s %s', T.Target(i), signStr);
    end

    if logical(T.BracketFound(i)) && ~isnan(T.C_instab_nM(i))
        C_instab(i) = sprintf('%.2f', double(T.C_instab_nM(i)));
        FoldShift(i) = sprintf('%.2f', double(T.FoldShift_vs_ref(i)));
        Trigger(i) = string(T.TriggerType(i));
    else
        C_instab(i) = sprintf('>%g or none found', double(T.TestedMax_nM(i)));
        FoldShift(i) = "--";
        Trigger(i) = "--";
    end

    if logical(T.BaseConverged(i)) && logical(T.BaseIsPhys(i))
        BaseStatus(i) = "usable";
    else
        BaseStatus(i) = "unusable";
    end
end

Tout = table(Perturbation, BaseStatus, C_instab, FoldShift, Trigger, ...
    'VariableNames', {'Perturbation','Baseline','C_instab_nM','Fold_vs_nominal','Trigger'});

disp(Tout)

csvOut = fullfile('results','oat_sensitivity_final_table.csv');
writetable(Tout, csvOut);
disp(csvOut)

texOut = fullfile('results','oat_sensitivity_final_table.tex');
fid = fopen(texOut, 'w');

fprintf(fid, '%% LaTeX-ready OAT sensitivity table\n');
fprintf(fid, '\\begin{table}[t]\n');
fprintf(fid, '\\caption{One-at-a-time sensitivity analysis in female endocardial T2DM. Each perturbation was applied to the T2DM layer, drug-free steady state was re-established, and $C_{\\mathrm{instab}}$ was re-evaluated.}\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{lllll}\n');
fprintf(fid, '\\hline\n');
fprintf(fid, 'Perturbation & Baseline & $C_{\\mathrm{instab}}$ (nM) & Fold vs nominal & Trigger \\\\\n');
fprintf(fid, '\\hline\n');

for i = 1:height(Tout)
    fprintf(fid, '%s & %s & %s & %s & %s \\\\\n', ...
        esc(Tout.Perturbation{i}), ...
        esc(Tout.Baseline{i}), ...
        esc(Tout.C_instab_nM{i}), ...
        esc(Tout.Fold_vs_nominal{i}), ...
        esc(Tout.Trigger{i}));
end

fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\label{tab:oat_sensitivity}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);

disp(texOut)

function s = esc(x)
    s = strrep(char(x), '_', '\_');
    s = strrep(s, '>', '$>$');
    s = strrep(s, '%', '\%%');
end