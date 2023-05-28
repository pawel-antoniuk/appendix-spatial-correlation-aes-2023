% Paweł Antoniuk 2023
% Bialystok University of Technology

%% Initialize
clearvars; clc;
addpath(genpath('SOFA'))
addpath(genpath('TwoEars-1.5'))
SOFAstart;
SOFAgetVersion()
% parpool(1)

%% Params
params.HRTFBaseDir = 'selHRTFs';
params.RecordingsBaseDir = 'recordings';
params.FinalResultsOutputDir = 'spatresults';
params.SpatOutputDir = params.FinalResultsOutputDir;
params.RecordingsExpectedFs = 48000;
params.RecordingSpatRange = [0 10];
params.RecordingFadeTime = [0.01 0.01];
params.RecordingLevelScale = 0.9;
params.NChannels = 2;
params.Elevations = 0;
params.IRmax = 3 * 512;
params.FadeDuration = 2*10^-3;
params.TargetTrackLoudness = -23; % db
params.MaxWidth = 0;
params.NRepetitions = 1;
params.AzimuthLocations = {
    [10 70]
    [10 70]
};
params.InverseAzimuthHRTFGroups = ["cipic"];
params.PostProcWin = hamming(1024 * 20);

%% Spatialize
geninputs(params)
spatialize(params)

%% Compute spatial correlations
poscParams.HRTFBaseDir = 'poscHRTFs';
poscParams.RecordingsExpectedFs = 48000;
poscParams.IRmax = 3 * 512;
poscParams.InverseAzimuthHRTFGroups = ["cipic"];
poscParams.FadeDuration = 2*10^-3;
poscParams.RMSThreshold = 0.1;
poscParams.AddNoiseRatio = 0.005;
poscParams.EnsembleWidthLimitThreshold = 0.7;
poscParams.POSCDistanceCorrection = 5;

HRTFs = loadHRTFs(poscParams);
[recordings, trackNames, realLocations] = loadRecordings(params);

targetPositions = linspace(-90, 90, 37*2);
targetPositions = [targetPositions' zeros(size(targetPositions, 2), 1)];
winWidth = size(params.PostProcWin, 1);
winOverlap = winWidth / 2;
winN = 1;
HRTFsN = length(HRTFs);
recordingsN = size(recordings, 1);
targetPositionsN = size(targetPositions, 1);
Cmat = zeros(HRTFsN, recordingsN, winN, targetPositionsN);
framesMask = false(HRTFsN, recordingsN, winN);

for iHRTF = 1:HRTFsN
    HRTF = HRTFs(iHRTF);
    HRTFSig = HRTF.SOFA.Data.IR ;
    HRTFPos = HRTF.Position;

    for iRecording = 1:recordingsN
        wholeSig = squeeze(recordings(iRecording, :, :));
        wholeSig = wholeSig + rand(size(wholeSig,1), 2) * poscParams.AddNoiseRatio;
        wholeSigRMSE = sqrt(mean(wholeSig .^ 2, 'all'));
        parfor iWindow = 1:winN
            winStart = winOverlap * (iWindow - 1) + 1;
            winStop = winStart + winWidth - 1;
            recSig = wholeSig;
            recSig(size(recSig, 1):winStop, :) = 0;
            recSig = recSig(winStart:winStop, :) .* params.PostProcWin;
            recSigRMS = sqrt(mean(recSig .^ 2, 'all'));

            if recSigRMS > poscParams.RMSThreshold * wholeSigRMSE
                framesMask(iHRTF, iRecording, iWindow) = true;
                for iPosition = 1:targetPositionsN
                    interpHRTF = interpolateHRTF( ...
                        HRTFSig, HRTFPos, targetPositions(iPosition, :), ...
                        Algorithm="VBAP");
                    hrtf = squeeze(interpHRTF)';
                    hrtf(size(hrtf,1)+1 : size(recSig,1), :) = 0;
                    hrtf = hrtf(1:size(recSig,1), :);
                    Cmat(iHRTF, iRecording, iWindow, iPosition) = ...
                        POSC(recSig, hrtf);
                end
            else 
                framesMask(iHRTF, iRecording, iWindow) = false;
            end
        end
    end
end

%% Correct POSCs
for iHRTF = 1:HRTFsN
    for iRecording = 2
        for iWindow = 1:winN
            if framesMask(iHRTF, iRecording, iWindow)
                for iPosition = 1:targetPositionsN
                    anglePos = targetPositions(iPosition);
                    correction = 1 + abs(anglePos) / 90 * poscParams.POSCDistanceCorrection;
                    Cmat(iHRTF, iRecording, iWindow, iPosition) = ...
                        Cmat(iHRTF, iRecording, iWindow, iPosition) * correction;
                end
            end
        end
    end
end


%% Calculate ensemble width limits
ensembleWidthLimits = zeros(HRTFsN, recordingsN, winN, 2);
for iHRTF = 1:HRTFsN
    for iRecording = 1:recordingsN
        for iWindow = 1:winN
            if framesMask(iHRTF, iRecording, iWindow)
                spatialCurve = squeeze(Cmat(iHRTF, iRecording, iWindow, :));
                threshold = max(spatialCurve) * poscParams.EnsembleWidthLimitThreshold;
                Istart = find(spatialCurve > threshold, 1, 'first');
                Iend = find(spatialCurve > threshold, 1, 'last');
                widthStart = Istart / targetPositionsN * 180 - 90;
                widhEnd = Iend / targetPositionsN * 180 - 90;
                ensembleWidthLimits(iHRTF, iRecording, iWindow, 1) = widthStart;
                ensembleWidthLimits(iHRTF, iRecording, iWindow, 2) = widhEnd;
            end
        end
    end
end

%% Calculate medians and ensemble widths
ensembleWidths = zeros(HRTFsN, recordingsN, 3);
for iHRTF = 1:HRTFsN
    for iRecording = 1:recordingsN
        mask = squeeze(framesMask(iHRTF, iRecording, :));
        selectedEnsembleWidthLimits = ...
            squeeze(ensembleWidthLimits(iHRTF, iRecording, mask, :));
        ensembleWidths(iHRTF, iRecording, 1) = ...
            median(ensembleWidthLimits(iHRTF, iRecording, :, 1));
        ensembleWidths(iHRTF, iRecording, 2) = ...
            median(ensembleWidthLimits(iHRTF, iRecording, :, 2));
        ensembleWidths(iHRTF, iRecording, 3) = ...
            ensembleWidths(iHRTF, iRecording, 2) ...
            - ensembleWidths(iHRTF, iRecording, 1);
    end
end

%% Calculate real ensemble widths

realEensembleWidthLimits = cellfun(@minmax, params.AzimuthLocations, ...
    'UniformOutput', false);
realEensembleWidths = cellfun(@(c) c(2) - c(1), realEensembleWidthLimits, ...
    'UniformOutput', false);

%% Plot
close all;
fig = figure();
fig.Position(3:4) = [600, 250];
tiledlayout(1, 2,"TileSpacing","tight")

for iRecording = 1:2
    recording = squeeze(recordings(iRecording, :, :));
    C = squeeze(sum(Cmat, 1));
    C = squeeze(C(iRecording, :, :));
    x = targetPositions(:, 1);

    nexttile
    plot(x, C')
    grid on
    xlabel("Azimuth \theta")
    xtickformat('%g°')
    if iRecording == 1
        ylabel('Spatial Correlation $C_\rho(\theta)$','Interpreter','latex')
    else
        ylabel('Spatial Correlation $\widetilde{C_\rho(\theta)}$','Interpreter','latex')
    end
    % ylim([-0.4, 0.9])
    xlim([0, 90])
    set(gca, 'XTick', -90:10:90)
    set(gca, 'YTick', -1:0.2:1)
    ax = gca;
    ax.TickDir = 'both';

    hold on
    plot([10 10], [-0.5 1], 'k--', LineWidth=1)
    plot([70 70], [-0.5 1], 'k--', LineWidth=1)
    text(11, 0.9, "white noise", HorizontalAlignment="left")
    text(69, 0.7, "white noise", HorizontalAlignment="right")
    hold off

end

exportgraphics(fig, "results/x1_whitenoise.png", Resolution=300);


%% Routines

function [recordings, trackNames, realLocations] = loadRecordings(params)
    filenames = dir(fullfile(params.SpatOutputDir, "*.wav"));

    firstRecording = audioread( ...
        fullfile(filenames(1).folder, ...
        filenames(1).name));
    recordings = zeros(length(filenames), size(firstRecording, 1), ...
        size(firstRecording, 2));
    recordings(1, :, :) = firstRecording;

    for iFilename = 2:length(filenames)
        recordings(iFilename, :, :) = audioread( ...
            fullfile(filenames(iFilename).folder, ...
            filenames(iFilename).name));
    end

    [~, filenames, ~] = fileparts({filenames.name});
    filenameParts = split(filenames, "_");

    if all(contains(filenameParts(1, :, 5), "loc"))
        locationParts = erase(filenameParts(1, :, 5), "loc");
        trackParts = filenameParts(1, :, 2);
    
        realLocations = cell(1, length(locationParts));
        trackNames = cell(1, length(trackParts));
    
        for iLocation = 1:length(locationParts)
            realLocations{iLocation} = str2double( ...
                split(locationParts{iLocation}, ','));
            trackNames{iLocation} = split(trackParts{iLocation}, ',');
        end
    
        [~, I] = sort(str2double(filenameParts(1, :, 1)));
        recordings = recordings(I, :, :);
        realLocations = realLocations(:, I);
    else
        trackNames = [];
        realLocations = str2double(erase(filenameParts(1, :, 4), "width"));
    end
end
