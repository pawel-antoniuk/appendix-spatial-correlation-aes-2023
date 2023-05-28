% Paweł Antoniuk 2023
% Bialystok University of Technology

%% Initialize
clearvars; clc;
addpath(genpath('SOFA'))
addpath(genpath('TwoEars-1.5'))
SOFAstart;
SOFAgetVersion()

%% Split train/test recordings
trainToAllRatio = 0.5;
recDirs = dir("input_recordings");
recDirs = recDirs([recDirs.isdir]);
recDirs = recDirs(3:end);
recordingNames = convertCharsToStrings({recDirs.name});
rng(10)

for ii = 1:3
    randRecordingNames = recordingNames(randperm(length(recordingNames)));
    trainI = 1:length(recordingNames) * trainToAllRatio;
    testI = length(recordingNames) * trainToAllRatio + 1 : length(recordingNames);
    trainTestRecordings(ii).Train = randRecordingNames(trainI);
    trainTestRecordings(ii).Test = randRecordingNames(testI);
end

fprintf("train recs: %d, test recs: %d\n", length(trainI), length(testI));

%% Params
params.HRTFBaseDir = 'selHRTFs';
params.FinalResultsOutputDir = 'frontNewMIT3';
params.SpatOutputDir = params.FinalResultsOutputDir;
params.RecordingsExpectedFs = 48000;
params.PostProcWin = hamming(1024);

%% Compute spatial correlations
poscParams.HRTFBaseDir = 'poscHRTFs';
poscParams.RecordingsExpectedFs = 48000;
poscParams.IRmax = 3 * 512;
poscParams.InverseAzimuthHRTFGroups = ["cipic"];
poscParams.FadeDuration = 2*10^-3;
poscParams.TargetPositions = linspace(-90, 90, 37);

%% Load data
HRTFs = loadHRTFs(poscParams);

%% Calculate train POSC
for iSplit = 1:length(trainTestRecordings)
    fprintf("Split %d / %d\n", iSplit, length(trainTestRecordings))

    [recordingsTrain, realWidthsTrain] = loadRecordings2( ...
        params, trainTestRecordings(iSplit).Train);
    [POSCsTrain, targetPositionsTrain] = generatePOSCs( ...
        params, poscParams, recordingsTrain, HRTFs);
    
    %% Predict train data
    spaceCorrection = linspace(0, 4, 10);
    spaceThreshold = linspace(0, 0.99, 10);
    results = cell(length(spaceCorrection), length(spaceThreshold));
    
    for iCorrection = 1:length(spaceCorrection)
        correction = spaceCorrection(iCorrection);
        for iThreshold = 1:length(spaceThreshold)
            tic
    
            threshold = spaceThreshold(iThreshold);
            predWidths = predictEnsembleWidths(...
                threshold, correction, POSCsTrain, ...
                targetPositionsTrain);
            predWidths = predWidths(:,:,3);
            mae = mean(abs(realWidthsTrain - predWidths));
            results{iCorrection, iThreshold}.RealWidths = realWidthsTrain;
            results{iCorrection, iThreshold}.PredWidths = predWidths;
            results{iCorrection, iThreshold}.Correction = correction;
            results{iCorrection, iThreshold}.Threshold = threshold;
            results{iCorrection, iThreshold}.Score = mae;
    
            fprintf("Pred %d / %d (%.2f s)\n", ...
                (iCorrection - 1) * length(spaceThreshold) + iThreshold, ...
                length(spaceCorrection) * length(spaceThreshold), ...
                toc)
        end
    end
    
    results = reshape(results, 1, []);
    results = [results{:}];
    
    %% Find best score
    bestResult = results(1);
    for iResult = 1:length(results)
        if results(iResult).Score < bestResult.Score
            bestResult = results(iResult);
        end
    end
    fprintf("Best: MAE = %.2f, threshold = %.2f, correction = %.2f\n", ...
        bestResult.Score, bestResult.Threshold, bestResult.Correction)
    
    %% Calculate test POSC
    [recordingsTest, realWidthsTest] = loadRecordings2( ...
        params, trainTestRecordings(iSplit).Test);
    [POSCsTest, targetPositionsTest] = generatePOSCs( ...
        params, poscParams, recordingsTest, HRTFs);
    
    %% Predict test data
    
    predWidths = predictEnsembleWidths(...
        bestResult.Threshold, bestResult.Correction, POSCsTest, ...
        targetPositionsTest);
    predWidths = predWidths(:, :, 3);

    finalResults(iSplit).RealTest = realWidthsTest;
    finalResults(iSplit).PredTest = predWidths;
    finalResults(iSplit).TrainRecordings = trainTestRecordings(iSplit).Train;
    finalResults(iSplit).TestRecordings = trainTestRecordings(iSplit).Test;
    finalResults(iSplit).TestMean = mean(abs(realWidthsTest - predWidths));
    finalResults(iSplit).TestStd = std(abs(realWidthsTest - predWidths));
    finalResults(iSplit).BestParams = bestResult;
    
    fprintf("Test MAE = %.2f\n", finalResults(iSplit).TestMean)
    fprintf("Test STD = %.2f\n", finalResults(iSplit).TestStd)
end

save(sprintf("workspace-dep-%s", datetime), '-v7.3')


function [POSCs, targetPositions] = generatePOSCs(params, poscParams, ...
    recordings, HRTFs)

targetPositions = poscParams.TargetPositions;
targetPositions = [targetPositions' zeros(size(targetPositions, 2), 1)];
winWidth = size(params.PostProcWin, 1);
winOverlap = winWidth / 2;
winN = floor(size(recordings, 2) / winOverlap);
HRTFsN = length(HRTFs);
recordingsN = size(recordings, 1);
targetPositionsN = size(targetPositions, 1);
POSCs = zeros(HRTFsN, recordingsN, winN, targetPositionsN);

for iHRTF = 1:HRTFsN
    HRTF = HRTFs(iHRTF);
    HRTFSig = HRTF.SOFA.Data.IR;
    HRTFPos = HRTF.Position;

    for iRecording = 1:recordingsN
        wholeSig = squeeze(recordings(iRecording, :, :));
        startTime = tic;
        parfor iWindow = 1:winN
            winStart = winOverlap * (iWindow - 1) + 1;
            winStop = winStart + winWidth - 1;
            recSig = wholeSig;
            recSig(size(recSig, 1):winStop, :) = 0;
            recSig = recSig(winStart:winStop, :) .* params.PostProcWin;
            for iPosition = 1:targetPositionsN
                interpHRTF = interpolateHRTF( ...
                    HRTFSig, HRTFPos, targetPositions(iPosition, :), ...
                    Algorithm="VBAP");
                hrtf = squeeze(interpHRTF)';
                hrtf(size(hrtf,1)+1 : size(recSig,1), :) = 0;
                hrtf = hrtf(1:size(recSig,1), :);
                POSCs(iHRTF, iRecording, iWindow, iPosition) = ...
                    POSC(recSig, hrtf);
            end
        end
            
        endTime = toc(startTime);

        fprintf("POSC %d / %d (%.2f s)\n", (iHRTF - 1) * HRTFsN + iRecording, ...
            HRTFsN * recordingsN, endTime)
    end
end
end


function ensembleWidths = predictEnsembleWidths(...
    thresholdRatio, correctionRatio, POSCs, targetPositions)
%% Correct POSCs
for iHRTF = 1:size(POSCs, 1)
    for iRecording = 1:size(POSCs, 2)
        for iWindow = 1:size(POSCs, 3)
            for iPosition = 1:size(POSCs, 4)
                anglePos = targetPositions(iPosition);
                correction = 1 + abs(anglePos) / 90 * correctionRatio;
                POSCs(iHRTF, iRecording, iWindow, iPosition) = ...
                    POSCs(iHRTF, iRecording, iWindow, iPosition) * correction;
            end
        end
    end
end

%% Calculate ensemble width limits
ensembleWidthLimits = zeros([size(POSCs, [1, 2, 3]), 2]);
for iHRTF = 1:size(POSCs, 1)
    for iRecording = 1:size(POSCs, 2)
        for iWindow = 1:size(POSCs, 3)
            spatialCurve = squeeze(POSCs(iHRTF, iRecording, iWindow, :));
            threshold = max(spatialCurve) * thresholdRatio;
            Istart = find(spatialCurve > threshold, 1, 'first');
            Iend = find(spatialCurve > threshold, 1, 'last');
            widthStart = Istart / size(POSCs, 4) * 180 - 90;
            widhEnd = Iend / size(POSCs, 4) * 180 - 90;
            ensembleWidthLimits(iHRTF, iRecording, iWindow, 1) = widthStart;
            ensembleWidthLimits(iHRTF, iRecording, iWindow, 2) = widhEnd;
        end
    end
end

%% Calculate medians and ensemble widths
ensembleWidths = zeros([size(POSCs, [1, 2]), 3]);
for iHRTF = 1:size(POSCs, 1)
    for iRecording = 1:size(POSCs, 2)
        ensembleWidths(iHRTF, iRecording, 1) = ...
            mean(ensembleWidthLimits(iHRTF, iRecording, :, 1));

        ensembleWidths(iHRTF, iRecording, 2) = ...
            mean(ensembleWidthLimits(iHRTF, iRecording, :, 2));

        ensembleWidths(iHRTF, iRecording, 3) = ...
            ensembleWidths(iHRTF, iRecording, 2) ...
            - ensembleWidths(iHRTF, iRecording, 1);

        ensembleWidths(iHRTF, iRecording, 3) = ...
            min(180, abs(ensembleWidths(iHRTF, iRecording, 3)));
    end
end
end

function nRecordings = getNumOfRecordings(params)
nRecordings = length(dir(fullfile(params.SpatOutputDir, "*.wav")));
end

function [recordings, realWidths, realLocs] = loadRecordings2(params, ...
    requestRecordingNames)
rng(10);
filenames = dir(fullfile(params.SpatOutputDir, "*.wav"));
filenames = filenames(randperm(length(filenames)));
filenameParts = split({filenames.name}, '_');
recordingNames = filenameParts(:, :, 1);
requestMask = cellfun(@(c) any(ismember(requestRecordingNames, c)), recordingNames);
filenames = filenames(requestMask);

firstRecording = audioread(fullfile(filenames(1).folder, ...
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
realWidths = erase(filenameParts(1, :, 4), 'width');
realWidths = str2double(convertCharsToStrings(realWidths));
realWidths = realWidths .* 2; % 0-90 -> 0-180
realLocs = erase(filenameParts(1, :, 6), 'azoffset');
realLocs = str2double(convertCharsToStrings(realLocs));
end
