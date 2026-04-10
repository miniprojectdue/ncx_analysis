function [flag, info] = detect_ead_last10beats(currents, bcl)
%DETECT_EAD_LAST10BEATS Detect EADs over the final 10 paced beats.
% EAD rule:
%   (i) after at least 50% repolarisation,
%   (ii) rebound >= 5 mV above preceding minimum,
%   (iii) max dV/dt > 10 V/s during the rebound.

    t = currents.time(:);
    V = currents.V(:);

    flag = false;
    info = struct('beatIndex', NaN, 'rebound_mV', NaN, 'max_dVdt', NaN);

    if isempty(t) || isempty(V)
        return;
    end

    tEnd = max(t);
    tStartWindow = max(min(t), tEnd - 10*bcl);

    beatStarts = tStartWindow:bcl:(tEnd - bcl + 1e-9);
    if isempty(beatStarts)
        beatStarts = max(min(t), tEnd - bcl);
    end

    for k = 1:numel(beatStarts)
        idx = t >= beatStarts(k) & t < beatStarts(k) + bcl;
        if nnz(idx) < 20
            continue;
        end

        tb = t(idx);
        Vb = V(idx);
        tb = tb - tb(1);

        % Ignore the very last sample if needed
        if numel(tb) < 10
            continue;
        end

        % Peak and repolarisation levels
        [Vmax, iPeak] = max(Vb);
        Vmin = min(Vb);
        ampV = Vmax - Vmin;

        if ampV < 40
            continue;
        end

        V50 = Vmax - 0.5 * ampV;

        % First time after peak when 50% repolarisation is reached
        i50rel = find(Vb(iPeak:end) <= V50, 1, 'first');
        if isempty(i50rel)
            continue;
        end
        i50 = iPeak + i50rel - 1;

        % Analyse only after 50% repolarisation
        tpost = tb(i50:end);
        Vpost = Vb(i50:end);

        if numel(tpost) < 5
            continue;
        end

        % Preceding minimum and rebound above it
        minSoFar = cummin(Vpost);
        rebound = Vpost - minSoFar;

        % dV/dt in V/s
        dVdt = gradient(Vpost, tpost) * 1000;  % mV/ms -> V/s

        reboundMask = rebound >= 5;
        if any(reboundMask)
            maxUpstroke = max(dVdt(reboundMask));
            maxRebound = max(rebound);

            if maxUpstroke > 10
                flag = true;
                info.beatIndex = k;
                info.rebound_mV = maxRebound;
                info.max_dVdt = maxUpstroke;
                return;
            end
        end
    end
end