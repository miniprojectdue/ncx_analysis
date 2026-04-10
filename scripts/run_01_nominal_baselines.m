setup_paths

celltypes = {'endo','epi','mid'};
sexes = {'male','female'};
diseases = {'control','t2dm'};

summaryRows = [];

for ic = 1:numel(celltypes)
    for is = 1:numel(sexes)
        for id = 1:numel(diseases)

            fprintf('Running %s | %s | %s\n', celltypes{ic}, sexes{is}, diseases{id});

            [param, X0] = build_base_param(celltypes{ic});
            param = apply_sex_layer(param, sexes{is});
            param = apply_t2dm_layer(param, diseases{id}, 'none');

            out = pace_to_convergence(X0, param, 500);
            fprintf('  converged=%d, isPhys=%d, beats=%d, RMP=%.3f, APD90=%.3f\n', ...
            out.converged, out.isPhys, out.nBeats, out.bm.RMP, out.bm.APD90);

            result = struct();
            result.celltype = celltypes{ic};
            result.sex = sexes{is};
            result.disease = diseases{id};
            result.param = param;
            result.bm = out.bm;
            result.currents = out.currents;
            result.X0_ss = out.X0;
            result.nBeats = out.nBeats;
            result.converged = out.converged;
            result.isPhys = out.isPhys;

            outFile = fullfile('results', sprintf('baseline_%s_%s_%s.mat', ...
                celltypes{ic}, sexes{is}, diseases{id}));
            save(outFile, 'result');
            summaryRows = [summaryRows; { ...
                celltypes{ic}, sexes{is}, diseases{id}, ...
                out.converged, out.isPhys, out.nBeats, ...
                out.bm.APD50, out.bm.APD90, out.bm.RMP, ...
                out.bm.CaT_amp, out.bm.CaD50, out.bm.CaD90}]; %#ok<AGROW>
            
        end
    end
end

summaryTable = cell2table(summaryRows, 'VariableNames', { ...
    'CellType','Sex','Disease','Converged','IsPhys','Beats', ...
    'APD50','APD90','RMP','CaT_amp','CaD50','CaD90'});

disp(summaryTable)
writetable(summaryTable, fullfile('results', 'baseline_summary.csv'));
save(fullfile('results', 'baseline_summary.mat'), 'summaryTable');