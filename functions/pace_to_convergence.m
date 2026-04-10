function out = pace_to_convergence(X0, param, maxBeats)
%PACE_TO_CONVERGENCE Pace in chunks until APD90 and CaT amplitude converge,
%and the final beat remains physiologically plausible.

    if nargin < 3
        maxBeats = 500;
    end

    options = [];
    chunk = 10;
    totalBeats = 0;

    apd90_hist = [];
    cat_hist = [];

    while totalBeats < maxBeats
        beats = min(chunk, maxBeats - totalBeats);

        [time, X] = modelRunner(X0, options, param, beats, 0);
        currents = extract_currents_wrapper(time, X, beats, param);

        X0 = X{end}(end, :);
        bm = extract_lastbeat_biomarkers(currents, param.bcl);

        apd90_hist(end+1) = bm.APD90; %#ok<AGROW>
        cat_hist(end+1) = bm.CaT_amp; %#ok<AGROW>
        totalBeats = totalBeats + beats;

        % Physiological validity check for the final beat
        isPhys = bm.RMP < -70 && bm.Vmax > 0 && bm.AP_amp > 80 && bm.APD90 < 0.9 * param.bcl;

        if numel(apd90_hist) >= 2
            dAPD = abs(apd90_hist(end) - apd90_hist(end-1));
            dCaT = abs(cat_hist(end) - cat_hist(end-1)) / max(abs(cat_hist(end-1)), eps);

            if dAPD < 1 && dCaT < 0.01 && isPhys
                out = struct();
                out.X0 = X0;
                out.currents = currents;
                out.bm = bm;
                out.nBeats = totalBeats;
                out.converged = true;
                out.isPhys = true;
                return;
            end
        end
    end

    % Final status if maxBeats reached
    isPhys = bm.RMP < -70 && bm.Vmax > 0 && bm.AP_amp > 80 && bm.APD90 < 0.9 * param.bcl;

    out = struct();
    out.X0 = X0;
    out.currents = currents;
    out.bm = bm;
    out.nBeats = totalBeats;
    out.converged = false;
    out.isPhys = isPhys;
end