function result = find_instability_threshold(X0_ss, baseParam, concGrid, maxDrugBeats)
%FIND_INSTABILITY_THRESHOLD Find first concentration producing repolarisation instability.
%
% Instability is defined as the first concentration at which any of:
%   (i) EAD is detected,
%   (ii) APD90 > 0.8*BCL,
%   (iii) physiological validity is lost.
%
% Output fields:
%   C_instab, beta_instab, bracketFound, triggerType, testedC, testedBeta, ...

    if nargin < 4
        maxDrugBeats = 200;
    end

    IC50 = 7;      % nM
    h = 0.65;

    testedC = [];
    testedBeta = [];
    testedInstab = [];
    testedTrigger = strings(0,1);
    testedConv = [];
    testedPhys = [];
    testedAPD90 = [];

    lowerC = NaN; lowerBeta = NaN;
    upperC = NaN; upperBeta = NaN;
    upperTrigger = "";

    prevNoInstab_C = NaN;
    prevNoInstab_B = NaN;

    % -------- coarse grid search --------
    for i = 1:numel(concGrid)
        C = concGrid(i);
        beta = (C^h) / (C^h + IC50^h);

        paramDrug = baseParam;
        paramDrug.IKr_Multiplier = getfield_def(paramDrug, 'IKr_Multiplier', 1) * (1 - beta);

        out = pace_to_convergence(X0_ss, paramDrug, maxDrugBeats);
        [eadFlag, ~] = detect_ead_last10beats(out.currents, paramDrug.bcl);

        apdFlag = ~isnan(out.bm.APD90) && out.bm.APD90 > 0.8 * paramDrug.bcl;
        physLossFlag = ~out.isPhys;

        instabFlag = eadFlag || apdFlag || physLossFlag;

        if eadFlag
            triggerType = "EAD";
        elseif apdFlag
            triggerType = "APD90>0.8BCL";
        elseif physLossFlag
            triggerType = "LossPhys";
        else
            triggerType = "None";
        end

        testedC(end+1,1) = C; %#ok<AGROW>
        testedBeta(end+1,1) = beta; %#ok<AGROW>
        testedInstab(end+1,1) = instabFlag; %#ok<AGROW>
        testedTrigger(end+1,1) = triggerType; %#ok<AGROW>
        testedConv(end+1,1) = out.converged; %#ok<AGROW>
        testedPhys(end+1,1) = out.isPhys; %#ok<AGROW>
        testedAPD90(end+1,1) = out.bm.APD90; %#ok<AGROW>

        if instabFlag
            upperC = C;
            upperBeta = beta;
            upperTrigger = triggerType;
            lowerC = prevNoInstab_C;
            lowerBeta = prevNoInstab_B;
            break;
        else
            prevNoInstab_C = C;
            prevNoInstab_B = beta;
        end
    end

    % No instability found on coarse grid
    if isnan(upperC)
        result = struct();
        result.C_instab = NaN;
        result.beta_instab = NaN;
        result.bracketFound = false;
        result.lowerC = NaN;
        result.upperC = NaN;
        result.triggerType = "";
        result.testedC = testedC;
        result.testedBeta = testedBeta;
        result.testedInstab = testedInstab;
        result.testedTrigger = testedTrigger;
        result.testedConv = testedConv;
        result.testedPhys = testedPhys;
        result.testedAPD90 = testedAPD90;
        return;
    end

    % If first point already unstable
    if isnan(lowerBeta)
        result = struct();
        result.C_instab = upperC;
        result.beta_instab = upperBeta;
        result.bracketFound = true;
        result.lowerC = NaN;
        result.upperC = upperC;
        result.triggerType = upperTrigger;
        result.testedC = testedC;
        result.testedBeta = testedBeta;
        result.testedInstab = testedInstab;
        result.testedTrigger = testedTrigger;
        result.testedConv = testedConv;
        result.testedPhys = testedPhys;
        result.testedAPD90 = testedAPD90;
        return;
    end

    % -------- beta bisection refinement --------
    maxIter = 10;
    relTol = 0.05;

    for it = 1:maxIter
        betaMid = 0.5 * (lowerBeta + upperBeta);
        Cmid = inverse_hill_conc(betaMid, IC50, h);

        paramDrug = baseParam;
        paramDrug.IKr_Multiplier = getfield_def(paramDrug, 'IKr_Multiplier', 1) * (1 - betaMid);

        out = pace_to_convergence(X0_ss, paramDrug, maxDrugBeats);
        [eadFlag, ~] = detect_ead_last10beats(out.currents, paramDrug.bcl);

        apdFlag = ~isnan(out.bm.APD90) && out.bm.APD90 > 0.8 * paramDrug.bcl;
        physLossFlag = ~out.isPhys;
        instabFlag = eadFlag || apdFlag || physLossFlag;

        if eadFlag
            triggerType = "EAD";
        elseif apdFlag
            triggerType = "APD90>0.8BCL";
        elseif physLossFlag
            triggerType = "LossPhys";
        else
            triggerType = "None";
        end

        testedC(end+1,1) = Cmid; %#ok<AGROW>
        testedBeta(end+1,1) = betaMid; %#ok<AGROW>
        testedInstab(end+1,1) = instabFlag; %#ok<AGROW>
        testedTrigger(end+1,1) = triggerType; %#ok<AGROW>
        testedConv(end+1,1) = out.converged; %#ok<AGROW>
        testedPhys(end+1,1) = out.isPhys; %#ok<AGROW>
        testedAPD90(end+1,1) = out.bm.APD90; %#ok<AGROW>

        if instabFlag
            upperBeta = betaMid;
            upperC = Cmid;
            upperTrigger = triggerType;
        else
            lowerBeta = betaMid;
            lowerC = Cmid;
        end

        if abs(upperC - lowerC) / max(upperC, eps) < relTol
            break;
        end
    end

    result = struct();
    result.C_instab = upperC;
    result.beta_instab = upperBeta;
    result.bracketFound = true;
    result.lowerC = lowerC;
    result.upperC = upperC;
    result.triggerType = upperTrigger;
    result.testedC = testedC;
    result.testedBeta = testedBeta;
    result.testedInstab = testedInstab;
    result.testedTrigger = testedTrigger;
    result.testedConv = testedConv;
    result.testedPhys = testedPhys;
    result.testedAPD90 = testedAPD90;
end

function val = getfield_def(s, fieldName, defaultVal)
    if isfield(s, fieldName)
        val = s.(fieldName);
    else
        val = defaultVal;
    end
end