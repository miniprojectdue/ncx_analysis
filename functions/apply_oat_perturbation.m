function [param, targetLabel] = apply_oat_perturbation(param, target, scale)
%APPLY_OAT_PERTURBATION Apply one-at-a-time perturbation to the T2DM layer.
%
% target: 'none', 'gkr', 'inaca', 'jup'
% scale : multiplicative factor, e.g. 0.8 or 1.2
%
% This function should be called AFTER apply_t2dm_layer(...,'t2dm',...)
% so that it perturbs the already-applied diseased multiplier.

    if nargin < 3
        scale = 1.0;
    end

    switch lower(target)
        case 'none'
            targetLabel = "Nominal";

        case 'gkr'
            if ~isfield(param, 'IKr_Multiplier')
                param.IKr_Multiplier = 1.0;
            end
            param.IKr_Multiplier = param.IKr_Multiplier * scale;
            targetLabel = "GKr";

        case 'inaca'
            if ~isfield(param, 'INaCa_Multiplier')
                param.INaCa_Multiplier = 1.0;
            end
            param.INaCa_Multiplier = param.INaCa_Multiplier * scale;
            targetLabel = "INaCa";

        case 'jup'
            if ~isfield(param, 'Jup_Multiplier')
                param.Jup_Multiplier = 1.0;
            end
            param.Jup_Multiplier = param.Jup_Multiplier * scale;
            targetLabel = "Jup";

        otherwise
            error('Unknown OAT target: %s', target);
    end
end