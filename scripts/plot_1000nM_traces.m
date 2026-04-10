setup_paths
clear functions
rehash

IC50 = 7;     % nM
h = 0.65;
C = 1000;     % nM
beta = (C^h) / (C^h + IC50^h);

% Conditions to inspect
cases = {
    'endo','male','control'
    'endo','male','t2dm'
    'endo','female','control'
    'endo','female','t2dm'
    'epi','female','control'
    'epi','male','t2dm'   % positive reference case
};

maxDrugBeats = 200;

for i = 1:size(cases,1)
    celltype = cases{i,1};
    sex      = cases{i,2};
    disease  = cases{i,3};

    fprintf('Plotting 1000 nM: %s | %s | %s\n', celltype, sex, disease);

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
    [eadFlag, eadInfo] = detect_ead_last10beats(out.currents, paramDrug.bcl);

    t = out.currents.time(:);
    V = out.currents.V(:);

    if isfield(out.currents, 'cai')
        Cai = out.currents.cai(:);
    elseif isfield(out.currents, 'Cai')
        Cai = out.currents.Cai(:);
    else
        Cai = [];
    end

    % Last 5 beats
    tEnd = max(t);
    idx5 = t >= (tEnd - 5*paramDrug.bcl);

    figure('Name', sprintf('%s_%s_%s_1000nM', celltype, sex, disease), ...
           'Color', 'w', 'Position', [100 100 1000 700]);

    subplot(3,1,1)
    plot(t(idx5), V(idx5), 'LineWidth', 1.2)
    xlabel('Time (ms)')
    ylabel('V_m (mV)')
    title(sprintf('%s | %s | %s | 1000 nM dofetilide | EAD=%d', ...
        celltype, sex, disease, eadFlag))
    grid on

    % Last beat only
    idx1 = t >= (tEnd - paramDrug.bcl);
    t1 = t(idx1);
    V1 = V(idx1);
    t1 = t1 - t1(1);

    subplot(3,1,2)
    plot(t1, V1, 'LineWidth', 1.2)
    xlabel('Time within last beat (ms)')
    ylabel('V_m (mV)')
    title(sprintf('Last beat | converged=%d | isPhys=%d | RMP=%.2f | APD90=%.2f', ...
        out.converged, out.isPhys, out.bm.RMP, out.bm.APD90))
    grid on

    % Calcium if available
    subplot(3,1,3)
    if ~isempty(Cai)
        plot(t(idx5), Cai(idx5), 'LineWidth', 1.2)
        xlabel('Time (ms)')
        ylabel('Cai')
        title('Last 5 beats: calcium')
        grid on
    else
        text(0.5, 0.5, 'No Cai field found', 'HorizontalAlignment', 'center')
        axis off
    end

    drawnow

    if eadFlag
        disp(eadInfo)
    end
end