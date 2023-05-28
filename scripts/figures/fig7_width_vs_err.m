preds = [finalResults(11,:).PredTest];
real = [finalResults(11,:).RealTest];
errs = preds - real;
win = 2;
x = linspace(0, 90, 1000);
resmea = zeros(1, length(x));
% resstd = zeros(1, length(x));

for ii = 1:length(x)
    mask = real >= x(ii) - win & real <= x(ii) + win;
    % resmea(ii) = mape(preds(mask), real(mask),"all");
    resmea(ii) = mean(abs(errs(mask)));
    % resstd(ii) = std(errs(real >= x(ii) - win & real <= x(ii) + win));
end

f = figure;
f.Position(3:4) = [400 300];

plot(x, resmea)

% hold on;
% plot(x, resstd)
% hold off;

xlim([0, 90])
% ylim([0, 30])
ylabel('Mean Absolute Error')
xlabel('Actual Ensemble Width \omega')

grid on
xtickformat('%g°')
ytickformat('%g°')
set(gca, 'XTick', -90:15:90)
% set(gca, 'YTick', -90:5:90)

exportgraphics(f, 'spat_actual_vs_err_dep.png', ...
    Resolution=300)