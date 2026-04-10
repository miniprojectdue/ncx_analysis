setup_paths
clear functions
rehash

Tthr = readtable(fullfile('results','threshold_summary.csv'), 'TextType', 'string');
Tinst = readtable(fullfile('results','instability_summary.csv'), 'TextType', 'string');
Tncx = readtable(fullfile('results','ncx_reset_followup.csv'), 'TextType', 'string');

% Compact final table rows
rowSpec = {
    'Endo','Male','Control','endo','male','control'
    'Endo','Male','T2DM','endo','male','t2dm'
    'Endo','Female','Control','endo','female','control'
    'Endo','Female','T2DM','endo','female','t2dm'
    'Epi','Male','Control','epi','male','control'
    'Epi','Male','T2DM','epi','male','t2dm'
};

n = size(rowSpec,1);

Cell   = strings(n,1);
Sex    = strings(n,1);
State  = strings(n,1);
C_EAD  = strings(n,1);
C_inst = strings(n,1);
Trigger = strings(n,1);
NCXreset = strings(n,1);

for i = 1:n
    Cell(i)  = rowSpec{i,1};
    Sex(i)   = rowSpec{i,2};
    State(i) = rowSpec{i,3};

    celltype = rowSpec{i,4};
    sex = rowSpec{i,5};
    disease = rowSpec{i,6};

    % --- strict EAD threshold ---
    idxThr = strcmpi(string(Tthr.CellType), celltype) & ...
             strcmpi(string(Tthr.Sex), sex) & ...
             strcmpi(string(Tthr.Disease), disease);

    if any(idxThr) && double(Tthr.BracketFound(find(idxThr,1))) == 1
        C_EAD(i) = sprintf('%.2f', double(Tthr.C_EAD_nM(find(idxThr,1))));
    else
        C_EAD(i) = "--";
    end

    % --- instability threshold ---
    idxInst = strcmpi(string(Tinst.CellType), celltype) & ...
              strcmpi(string(Tinst.Sex), sex) & ...
              strcmpi(string(Tinst.Disease), disease);

    if any(idxInst) && double(Tinst.BracketFound(find(idxInst,1))) == 1
        C_inst(i) = sprintf('%.2f', double(Tinst.C_instab_nM(find(idxInst,1))));
        trig = string(Tinst.TriggerType(find(idxInst,1)));
        if strlength(trig) == 0
            Trigger(i) = "--";
        else
            Trigger(i) = trig;
        end
    else
        C_inst(i) = ">1000";
        Trigger(i) = "--";
    end

    % --- NCX reset column only for endocardial T2DM rows ---
    if strcmpi(celltype,'endo') && strcmpi(disease,'t2dm')
        idxNCX = strcmpi(string(Tncx.CellType), celltype) & ...
                 strcmpi(string(Tncx.Sex), sex) & ...
                 strcmpi(string(Tncx.Variant), 'T2DM_NCXreset');

        if any(idxNCX) && double(Tncx.BracketFound(find(idxNCX,1))) == 1
            NCXreset(i) = sprintf('%.2f', double(Tncx.C_instab_nM(find(idxNCX,1))));
        else
            NCXreset(i) = ">1000";
        end
    else
        NCXreset(i) = "--";
    end
end

finalTable = table(Cell, Sex, State, C_EAD, C_inst, Trigger, NCXreset, ...
    'VariableNames', {'Cell','Sex','State','C_EAD_nM','C_instab_nM','Trigger','NCXreset_C_instab_nM'});

csvOut = fullfile('results','final_results_table.csv');
writetable(finalTable, csvOut)
disp(csvOut)

% --- write LaTeX table fragment ---
texOut = fullfile('results','final_results_table.tex');
fid = fopen(texOut, 'w');

fprintf(fid, '%% LaTeX-ready table fragment\n');
fprintf(fid, '\\begin{table}[t]\n');
fprintf(fid, '\\caption{Summary of strict EAD and broader instability thresholds. Concentrations are in nM. For endocardial T2DM rows, the final column reports the instability threshold after resetting the T2DM NCX multiplier to control level.}\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{lllllll}\n');
fprintf(fid, '\\hline\n');
fprintf(fid, 'Cell & Sex & State & $C_{\\mathrm{EAD}}$ & $C_{\\mathrm{instab}}$ & Trigger & NCX-reset \\\\\n');
fprintf(fid, '\\hline\n');

for i = 1:height(finalTable)
    fprintf(fid, '%s & %s & %s & %s & %s & %s & %s \\\\\n', ...
        esc(finalTable.Cell{i}), ...
        esc(finalTable.Sex{i}), ...
        esc(finalTable.State{i}), ...
        esc(finalTable.C_EAD_nM{i}), ...
        esc(finalTable.C_instab_nM{i}), ...
        esc(finalTable.Trigger{i}), ...
        esc(finalTable.NCXreset_C_instab_nM{i}));
end

fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\label{tab:threshold_summary}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);

disp(texOut)

function s = esc(x)
    s = strrep(char(x), '_', '\_');
    s = strrep(s, '>', '$>$');
end