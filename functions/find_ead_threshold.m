function result = find_ead_threshold(X0_ss, baseParam, concGrid, maxDrugBeats)
%FIND_EAD_THRESHOLD Bracket and refine dofetilide EAD threshold.
%
% Inputs:
%   X0_ss        - drug-free steady-state vector for this condition
%   baseParam    - parameter struct for this condition
%   concGrid     - coarse nM grid, e.g. [0.1 0.3 1 3 10 30 100 300 1000]
%   maxDrugBeats - pacing cap after drug application, e.g. 200
%
% Output:
%   result struct with C_EAD, beta_EAD, bracket, tested concentrations, etc.

    if nargin < 4
        maxDrugBeats = 200;
    end

    IC50 = 7;      % nM
    h = 0.65;

    testedC = [];
    testedBeta = [];
    testedEAD = [];
    testedConv = [];
    testedPhys = [];

    lowerC = NaN; lowerBeta = NaN;
    upperC = NaN; upperBeta = NaN;

    % -------- coarse grid search --------
    prevNoEAD_C = NaN;
    prevNoEAD_B = NaN;

    for i = 1:numel(concGrid)
        C = concGrid(i);
        beta = (C^h) / (C^h + IC50^h);

        paramDrug = baseParam;
        paramDrug.IKr_Multiplier = getfield_def(paramDrug, 'IKr_Multiplier', 1) * (1 - beta);

        out = pace_to_convergence(X0_ss, paramDrug, maxDrugBeats);
        [eadFlag, eadInfo] = detect_ead_last10beats(out.currents, paramDrug.bcl); %#ok<NASGU>

        testedC(end+1,1) = C; %#ok<AGROW>
        testedBeta(end+1,1) = beta; %#ok<AGROW>
        testedEAD(end+1,1) = eadFlag; %#ok<AGROW>
        testedConv(end+1,1) = out.converged; %#ok<AGROW>
        testedPhys(end+1,1) = out.isPhys; %#ok<AGROW>

        if eadFlag
            upperC = C;
            upperBeta = beta;
            lowerC = prevNoEAD_C;
            lowerBeta = prevNoEAD_B;
            break;
        else
            prevNoEAD_C = C;
            prevNoEAD_B = beta;
        end
    end

    % No EAD found on coarse grid
    if isnan(upperC)
        result = struct();
        result.C_EAD = NaN;
        result.beta_EAD = NaN;
        result.bracketFound = false;
        result.lowerC = NaN;
        result.upperC = NaN;
        result.testedC = testedC;
        result.testedBeta = testedBeta;
        result.testedEAD = testedEAD;
        result.testedConv = testedConv;
        result.testedPhys = testedPhys;
        return;
    end

    % If first grid point already had EAD, stop there
    if isnan(lowerBeta)
        result = struct();
        result.C_EAD = upperC;
        result.beta_EAD = upperBeta;
        result.bracketFound = true;
        result.lowerC = NaN;
        result.upperC = upperC;
        result.testedC = testedC;
        result.testedBeta = testedBeta;
        result.testedEAD = testedEAD;
        result.testedConv = testedConv;
        result.testedPhys = testedPhys;
        return;
    end

    % -------- beta bisection refinement --------
    maxIter = 10;
    relTol = 0.05;  % 5%

    for it = 1:maxIter
        betaMid = 0.5 * (lowerBeta + upperBeta);
        Cmid = inverse_hill_conc(betaMid, IC50, h);

        paramDrug = baseParam;
        paramDrug.IKr_Multiplier = getfield_def(paramDrug, 'IKr_Multiplier', 1) * (1 - betaMid);

        out = pace_to_convergence(X0_ss, paramDrug, maxDrugBeats);
        [eadFlag, eadInfo] = detect_ead_last10beats(out.currents, paramDrug.bcl); %#ok<NASGU>

        testedC(end+1,1) = Cmid; %#ok<AGROW>
        testedBeta(end+1,1) = betaMid; %#ok<AGROW>
        testedEAD(end+1,1) = eadFlag; %#ok<AGROW>
        testedConv(end+1,1) = out.converged; %#ok<AGROW>
        testedPhys(end+1,1) = out.isPhys; %#ok<AGROW>

        if eadFlag
            upperBeta = betaMid;
            upperC = Cmid;
        else
            lowerBeta = betaMid;
            lowerC = Cmid;
        end

        if abs(upperC - lowerC) / max(upperC, eps) < relTol
            break;
        end
    end

    result = struct();
    result.C_EAD = upperC;
    result.beta_EAD = upperBeta;
    result.bracketFound = true;
    result.lowerC = lowerC;
    result.upperC = upperC;
    result.testedC = testedC;
    result.testedBeta = testedBeta;
    result.testedEAD = testedEAD;
    result.testedConv = testedConv;
    result.testedPhys = testedPhys;
end

function val = getfield_def(s, fieldName, defaultVal)
    if isfield(s, fieldName)
        val = s.(fieldName);
    else
        val = defaultVal;
    end
end