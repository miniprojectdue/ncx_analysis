function bm = extract_lastbeat_biomarkers(currents, bcl)
%EXTRACT_LASTBEAT_BIOMARKERS Extract simple last-beat AP/Ca biomarkers.

    t = currents.time(:);
    V = currents.V(:);

    if isfield(currents, 'cai')
        Cai = currents.cai(:);
    elseif isfield(currents, 'Cai')
        Cai = currents.Cai(:);
    else
        error('Could not find cai in currents structure.');
    end

    tEnd = max(t);
    idx = t >= (tEnd - bcl);

    tb = t(idx);
    Vb = V(idx);
    Caib = Cai(idx);

    tb = tb - tb(1);

    Vmax = max(Vb);
    Vmin = min(Vb);
    ampV = Vmax - Vmin;

    level90 = Vmax - 0.90 * ampV;
    level50 = Vmax - 0.50 * ampV;

    iPeak = find(Vb == Vmax, 1, 'first');
    postV = Vb(iPeak:end);
    postt = tb(iPeak:end);

    i90 = find(postV <= level90, 1, 'first');
    i50 = find(postV <= level50, 1, 'first');

    if isempty(i90)
        APD90 = NaN;
    else
        APD90 = postt(i90) - tb(iPeak);
    end

    if isempty(i50)
        APD50 = NaN;
    else
        APD50 = postt(i50) - tb(iPeak);
    end

    Ca_min = min(Caib);
    Ca_max = max(Caib);
    Ca_amp = Ca_max - Ca_min;

    levelCa50 = Ca_min + 0.50 * Ca_amp;
    levelCa90 = Ca_min + 0.10 * Ca_amp;

    iCaPeak = find(Caib == Ca_max, 1, 'first');
    postCa = Caib(iCaPeak:end);
    postCa_t = tb(iCaPeak:end);

    iCa50 = find(postCa <= levelCa50, 1, 'first');
    iCa90 = find(postCa <= levelCa90, 1, 'first');

    if isempty(iCa50)
        CaD50 = NaN;
    else
        CaD50 = postCa_t(iCa50) - tb(iCaPeak);
    end

    if isempty(iCa90)
        CaD90 = NaN;
    else
        CaD90 = postCa_t(iCa90) - tb(iCaPeak);
    end
    AP_amp = Vmax - Vmin;
    bm = struct();
    bm.APD50 = APD50;
    bm.APD90 = APD90;
    bm.AP_amp = AP_amp;
    bm.RMP = Vb(1);
    bm.Vmax = Vmax;
    bm.CaT_amp = Ca_amp;
    bm.Ca_peak = Ca_max;
    bm.Ca_diastolic = Caib(1);
    bm.CaD50 = CaD50;
    bm.CaD90 = CaD90;
end