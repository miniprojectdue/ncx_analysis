function [param, X0] = build_base_param(celltype)
%BUILD_BASE_PARAM Set base pacing/model settings and starting state.

    param = struct();
    param.bcl = 1000;
    param.verbose = false;
    param.stimAmp = -53;
    param.stimDur = 1;
    param.cellType = [];   % 0 endo, 1 epi, 2 mid

    switch lower(celltype)
        case 'endo'
            param.cellType = 0;
            X0 = getStartingState('m_endo');
        case 'epi'
            param.cellType = 1;
            X0 = getStartingState('m_epi');
        case 'mid'
            param.cellType = 2;
            X0 = getStartingState('m_mid');
        otherwise
            error('Unknown cell type: %s', celltype);
    end
end