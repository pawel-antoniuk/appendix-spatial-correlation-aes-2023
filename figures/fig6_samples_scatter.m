T = readtable('actual_vs_predicted.csv');

f = figure;
f.Position(3:4) = [400 300];

scatter(T.predicted, T.actual, 20, 'filled', ...
    MarkerFaceAlpha=0.1, ...
    MarkerEdgeAlpha=0.1)
xlim([0, 180])
ylim([0, 180])
ylabel('Actual Ensemble Width')
xlabel('Predicted Ensemble Width')

hold on

plot([0, 180], [0, 180], LineWidth=2)
% 
% p = polyfit(T.predicted, T.actual, 5);
% px = linspace(0, 90);
% py = polyval(p, px);
% plot(px, py)

hold off

legend('Samples', 'One to one line', ...
    Location='southeast')
grid on

exportgraphics(f, 'actual_vs_predicted.png', ...
    Resolution=300)