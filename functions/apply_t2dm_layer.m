function param = apply_t2dm_layer(param, disease, ical_branch)
%APPLY_T2DM_LAYER Apply control or T2DM multipliers.
%
% disease: 'control' or 't2dm'
% ical_branch: 'none', 'up', or 'down'

    if nargin < 3
        ical_branch = 'none';
    end

    if strcmpi(disease, 'control')
        return;
    end

    param.IKr_Multiplier   = getfield_def(param, 'IKr_Multiplier',   1) * 0.35;
    param.IKs_Multiplier   = getfield_def(param, 'IKs_Multiplier',   1) * 0.95;
    param.INaCa_Multiplier = getfield_def(param, 'INaCa_Multiplier', 1) * 3.44;
    param.Jup_Multiplier   = getfield_def(param, 'Jup_Multiplier',   1) * 0.70;

    switch lower(ical_branch)
        case 'none'
        case 'up'
            param.ICaL_Multiplier = getfield_def(param, 'ICaL_Multiplier', 1) * 1.14;
        case 'down'
            param.ICaL_Multiplier = getfield_def(param, 'ICaL_Multiplier', 1) * 0.62;
        otherwise
            error('Unknown ICaL branch: %s', ical_branch);
    end
end

function val = getfield_def(s, fieldName, defaultVal)
    if isfield(s, fieldName)
        val = s.(fieldName);
    else
        val = defaultVal;
    end
end