setup_paths
clear functions
rehash

Tin = readtable(fullfile('results','oat_sensitivity_summary.csv'), 'TextType', 'string');

idx = strcmpi(Tin.CellType, 'endo') & strcmpi(Tin.Sex, 'female');
T = Tin(idx,:);

orderLabels = [
    "Nominal_x1.00"
    "GKr_x0.80"
    "GKr_x1.20"
    "INaCa_x0.80"
    "INaCa_x1.20"
    "Jup_x0.80"
    "Jup_x1.20"
];

rowKey = strings(height(T),1);
for i = 1:height(T)
    rowKey(i) = sprintf('%s_x%.2f', string(T.Target(i)), double(T.Scale(i)));
end

[~, loc] = ismember(orderLabels, rowKey);
loc = loc(loc > 0);
T = T(loc,:);

labels = strings(height(T),1);
yvals = nan(height(T),1);

for i = 1:height(T)
    if strcmpi(T.Target(i), "Nominal")
        labels(i) = "Nominal";
    else
        signStr = "+20%";
        if double(T.Scale(i)) < 1
            signStr = "-20%";
        end
        labels(i) = sprintf('%s %s', T.Target(i), signStr);
    end

    if logical(T.BracketFound(i)) && ~isnan(T.C_instab_nM(i))
        yvals(i) = double(T.FoldShift_vs_ref(i));
    end
end

fig = figure('Color','w','Position',[100 100 900 420]);
ax = axes(fig);
hold(ax, 'on')

x = 1:numel(yvals);
plot(ax, x, yvals, 'o', 'MarkerSize', 8, 'LineWidth', 1.2)
plot(ax, x, yvals, '-', 'LineWidth', 1.0)

yline(ax, 1.0, '--', 'Nominal', 'LineWidth', 1.0);

set(ax, 'XTick', x, 'XTickLabel', labels)
xtickangle(ax, 25)
ylabel(ax, 'C_{instab} / nominal C_{instab}')
title(ax, 'Female endocardial T2DM: OAT sensitivity of instability threshold')
grid(ax, 'on')
box(ax, 'off')

outPng = fullfile('results','oat_sensitivity_foldshift.png');
outPdf = fullfile('results','oat_sensitivity_foldshift.pdf');

exportgraphics(fig, outPng, 'Resolution', 300)
exportgraphics(fig, outPdf, 'ContentType', 'vector')

disp(outPng)
disp(outPdf)