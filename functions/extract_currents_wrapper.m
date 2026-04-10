function currents = extract_currents_wrapper(time, X, beats, param)
%EXTRACT_CURRENTS_WRAPPER Repo-specific helper wrapper.
    currents = getCurrentsStructure(time, X, beats, param, 0);
end