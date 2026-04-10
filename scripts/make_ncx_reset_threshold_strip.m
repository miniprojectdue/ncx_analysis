setup_paths
clear functions
rehash

fig = figure('Color','w','Position',[100 100 900 420]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

xmax = 1000;
xmin = 1;

ypos  = [1 2 3];
ylabs = {'T2DM + NCX reset','T2DM','Control'};

% -------------------------
% Panel 1: Male endocardium
% -------------------------
nexttile
hold on

xvals = [xmax, 90.96, xmax];
isCensored = [true, false, true];

for i = 1:numel(ypos)
    y = ypos(i);

    plot([xmin xmax], [y y], ':', 'Color', [0.82 0.82 0.82], 'LineWidth', 1);

    if isCensored(i)
        plot(xmax, y, 'o', ...
            'MarkerSize', 7, ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.2);
        quiver(xmax*0.93, y, xmax*0.10, 0, 0, ...
            'Color', 'k', 'LineWidth', 1.2, 'MaxHeadSize', 0.8);
        text(xmax/1.16, y+0.14, '>1000', ...
            'FontSize', 8, 'HorizontalAlignment', 'center');
    else
        plot(xvals(i), y, 'o', ...
            'MarkerSize', 8, ...
            'MarkerFaceColor', 'k', ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.2);
        text(xvals(i), y+0.14, sprintf('%.2f', xvals(i)), ...
            'FontSize', 9, 'HorizontalAlignment', 'center');
    end
end

set(gca, 'XScale', 'log')
xlim([xmin xmax*1.18])
ylim([0.5 3.5])
set(gca, 'YTick', ypos, 'YTickLabel', ylabs)
ylabel('Male endo')
grid on
box off

% ---------------------------
% Panel 2: Female endocardium
% ---------------------------
nexttile
hold on

xvals = [xmax, 4.25, xmax];
isCensored = [true, false, true];

for i = 1:numel(ypos)
    y = ypos(i);

    plot([xmin xmax], [y y], ':', 'Color', [0.82 0.82 0.82], 'LineWidth', 1);

    if isCensored(i)
        plot(xmax, y, 'o', ...
            'MarkerSize', 7, ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.2);
        quiver(xmax*0.93, y, xmax*0.10, 0, 0, ...
            'Color', 'k', 'LineWidth', 1.2, 'MaxHeadSize', 0.8);
        text(xmax/1.16, y+0.14, '>1000', ...
            'FontSize', 8, 'HorizontalAlignment', 'center');
    else
        plot(xvals(i), y, 'o', ...
            'MarkerSize', 8, ...
            'MarkerFaceColor', 'k', ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.2);
        text(xvals(i), y+0.14, sprintf('%.2f', xvals(i)), ...
            'FontSize', 9, 'HorizontalAlignment', 'center');
    end
end

set(gca, 'XScale', 'log')
xlim([xmin xmax*1.18])
ylim([0.5 3.5])
set(gca, 'YTick', ypos, 'YTickLabel', ylabs)
ylabel('Female endo')
xlabel('Dofetilide concentration, C_{instab} (nM)')
grid on
box off

sgtitle('Endocardial instability thresholds with and without NCX reset')

outPng = fullfile('results','ncx_reset_threshold_strip.png');
outPdf = fullfile('results','ncx_reset_threshold_strip.pdf');

exportgraphics(fig, outPng, 'Resolution', 300)
exportgraphics(fig, outPdf, 'ContentType', 'vector')

disp(outPng)
disp(outPdf)