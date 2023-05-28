load('workspace-indep2-30-Mar-2023 07_25_34.mat')
location = [finalResults(:,:).RealLocsTest]';
width = [finalResults(:,:).RealTest]';
err = ([finalResults(:,:).RealTest] - [finalResults(:,:).PredTest])';

locationSpace = linspace(-45, 45, 200);
widthSpace = linspace(0, 90, 200);
winWidth = 15;
map = zeros(length(locationSpace), length(widthSpace));

for iLocation = 1:length(locationSpace)
    locationStart = locationSpace(iLocation) - winWidth / 2;
    locationEnd = locationSpace(iLocation) + winWidth / 2;

    for iWidth = 1:length(widthSpace)
        widthStart = widthSpace(iWidth) - winWidth / 2;
        widthEnd = widthSpace(iWidth) + winWidth / 2;

        errors = err(location >= locationStart & location < locationEnd ...
            & width >= widthStart & width < widthEnd, :);

%         if sum(errors) > 0
            mae = mean(abs(errors));
            map(iLocation, iWidth) = mae;
%         end
    end
end
f = figure;
f.Position(3:4) = [400 300];

imagesc(locationSpace, widthSpace, map')
set(gca,'YDir','normal')
set(gca, 'XTick', -45:15:45, 'XTickLabel', -45:15:45);
set(gca, 'YTick', -90:15:90, 'YTickLabel', -90:15:90);
xtickformat('%g°')
ytickformat('%g°')
ylabel('Ensemble Width \omega')
xlabel('Ensemble Location \phi')
set(gca,'TickDir', 'both')
grid on
% set(gca, 'LineWidth', 1)
set(gca, 'GridAlpha', 1)
colormap turbo
ylim([0 90])
xlim([-45 45])

colormap jet

c = colorbar;
c.Label.String = 'Mean Absolute Error';
c.Ruler.TickLabelFormat = '%g°';

% traingle
hold on
p = polyshape([-45 0 45],[60 15 60]);
pg = plot(p);
pg.FaceAlpha = 0;
pg.LineStyle = '--';
% pg.EdgeColor = [1 1 1];
pg.LineWidth = 1;
hold off

exportgraphics(f, 'spat_mae_loc_width.png', ...
    Resolution=300)


