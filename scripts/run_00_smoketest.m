setup_paths

clear param
param.bcl = 1000;
param.model = @model_ToRORd_Land_Male;
param.verbose = true;
param.cellType = 0;   % 0 = endo, 1 = epi, 2 = mid

options = [];
beats = 20;
ignoreFirst = beats - 1;

X0 = getStartingState('m_endo');

[time, X] = modelRunner(X0, options, param, beats, ignoreFirst);
currents = getCurrentsStructure(time, X, beats, param, 0);

figure;
plot(currents.time, currents.V, 'LineWidth', 1.2)
xlabel('Time (ms)')
ylabel('V_m (mV)')
title('Smoke test: male endocardial')
grid on