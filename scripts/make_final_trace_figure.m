setup_paths
clear functions
rehash

IC50 = 7;     % nM
h = 0.65;
C = 1000;     % nM
beta = (C^h) / (C^h + IC50^h);
maxDrugBeats = 200;

% Final figure panel order
cases = {
    'endo','male','control'
    'endo','male','t2dm'
    'endo','female','control'
    'endo','female','t2dm'
    'epi','male','control'
    'epi','male','t2dm'
};

fig = figure('Color','w','Position',[100 100 1200 700]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

for i = 1:size(cases,1)
    celltype = cases{i,1};
    sex      = cases{i,2};
    disease  = cases{i,3};

    inFile = fullfile('results', sprintf('baseline_%s_%s_%s.mat', celltype, sex, disease));
    S = load(inFile, 'result');

    X0_ss = S.result.X0_ss;
    baseParam = S.result.param;

    paramDrug = baseParam;
    if isfield(paramDrug, 'IKr_Multiplier')
        paramDrug.IKr_Multiplier = paramDrug.IKr_Multiplier * (1 - beta);
    else
        paramDrug.IKr_Multiplier = (1 - beta);
    end

    out = pace_to_convergence(X0_ss, paramDrug, maxDrugBeats);
    [eadFlag, ~] = detect_ead_last10beats(out.currents, paramDrug.bcl);

    t = out.currents.time(:);
    V = out.currents.V(:);

    idx = t >= (max(t) - paramDrug.bcl);
    tb = t(idx);
    Vb = V(idx);
    tb = tb - tb(1);

    nexttile
    plot(tb, Vb, 'LineWidth', 1.2)
    xlabel('Time within beat (ms)')
    ylabel('V_m (mV)')
    title(sprintf('%s | %s | %s', celltype, sex, disease), 'Interpreter', 'none')
    subtitle(sprintf('EAD=%d, APD90=%.2f ms', eadFlag, out.bm.APD90), 'Interpreter', 'none')
    xlim([0 1000])
    ylim([-95 40])
    grid on
end

sgtitle('Representative last-beat traces at 1000 nM dofetilide')

outPng = fullfile('results','final_trace_figure_1000nM.png');
outPdf = fullfile('results','final_trace_figure_1000nM.pdf');

exportgraphics(fig, outPng, 'Resolution', 300)
exportgraphics(fig, outPdf, 'ContentType', 'vector')
disp(outPng)
disp(outPdf)