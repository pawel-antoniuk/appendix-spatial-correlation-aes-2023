% Paweł Antoniuk 2023
% Bialystok University of Technology

%% Initialize
clearvars; close all; clc;
addpath('SOFA/API_MO');
SOFAstart;
SOFAgetVersion()
% parpool(8)

%% Params
params.HRTFBaseDir = 'smallHRTFs';
params.RecordingsBaseDir = 'input_recordings';
params.FinalResultsOutputDir = 'rec_front_small_uniform';
params.SpatOutputDir = [params.FinalResultsOutputDir filesep 'spat'];
params.RecordingsExpectedFs = 48000;
params.RecordingLoadRange = [2.5 inf];
params.RecordingSpatRange = [0.5 7];
params.RecordingFadeTime = [0.01 0.01];
params.RecordingLevelScale = 0.9;
params.NChannels = 2;
params.Elevations = [0];
params.IRmax = 3 * 512;
params.FadeDuration = 2*10^-3;
params.TargetTrackLoudness = -23; % db
params.MaxWidth = 45; % from ensemble center
params.NRepetitions = 8;
params.AzimuthOffsetRange = [-45 45];
params.InverseAzimuthHRTFGroups = ["cipic"];

%% Load HRTFs and get audio filenames
HRTFs = loadHRTFs(params);
audioFilenames = dirWithoutDots(params.RecordingsBaseDir);

%% Spatialize songs
nAudioFilenames = length(audioFilenames);
allSpatMetaresults = cell(nAudioFilenames, params.NRepetitions);

start = tic;
for iRep = 1:params.NRepetitions
    % parfor
    parfor iAudioFilename = 1:nAudioFilenames       
        %% Load audio tracks and spatialize them 
        audioFilename = audioFilenames(iAudioFilename);
        [tracks,trackNames] = loadAudioTracks(audioFilename, params);
        tracks = normalizeAudioTracks(tracks, params);
        [spatResults,spatMetaresults] = spatializeSong(HRTFs, tracks, ...
            trackNames, audioFilename.name, params);
        allSpatMetaresults{iAudioFilename,iRep} = spatMetaresults;
        
        %% Postprocess results
        spatResults = posprocessSpatResults(spatResults, params);
        
        %% Save results
        saveSpatResults(spatResults, spatMetaresults, HRTFs, ...
            audioFilename, iRep, params);
        
        fprintf("Progress  [audio: %d/%d] (%s)\n", ...
            iAudioFilename, nAudioFilenames, audioFilename.name);
    end
end
toc(start)

%% Plot scenes summary
% spatMetaresults = cell2mat(reshape(allSpatMetaresults, 1, 1, 1, []));
% plotAudioScene(HRTFs, spatMetaresults, params);

%% Save workspace
save(fullfile(params.FinalResultsOutputDir, 'workspace'), '-v7.3');
saveHRTFMetadata(HRTFs, params);

%% --- ROUTINES ---
% Load HRTFs routine
function HRTFs = loadHRTFs(params)
    HRTFFilenames = dir(fullfile(params.HRTFBaseDir, '*', '*.sofa'));
    
    % HRTF struct definition
    HRTFs = struct('Id', [], ...
        'Name', [], ...
        'Folder', [], ...
        'HRTFGroup', [], ...
        'SOFA', [], ...
        'Position', [], ...
        'Distance', []);
    HRTFGroupData = containers.Map;
    
    for iHRTF = 1:length(HRTFFilenames)
        filename = HRTFFilenames(iHRTF);
        fullFilename = fullfile(filename.folder, filename.name);
        
        HRTFs(iHRTF) = loadHRTF(iHRTF, fullFilename, params);   
        
        if HRTFs(iHRTF).SOFA.Data.SamplingRate ~= params.RecordingsExpectedFs
            [loadStatus,HRTFs(iHRTF)] = tryLoadResampledHRTF(iHRTF, ...
                HRTFs(iHRTF), params);
            if ~loadStatus
                resampleAndSave(HRTFs(iHRTF), params);
                [loadStatus,HRTFs(iHRTF)] = tryLoadResampledHRTF(iHRTF, ...
                    HRTFs(iHRTF), params);
                
                if ~loadStatus
                    error('Cannot find previously resampled HRTF');
                end
            end
        end

        ir =  HRTFs(iHRTF).SOFA.Data.IR;
        if size(ir, 3) > params.IRmax
            Nfadeout = params.FadeDuration*params.RecordingsExpectedFs;
            fade = [repelem(1,params.IRmax-Nfadeout), ...
                (cos(linspace(0,pi,Nfadeout))+1)/2];
            fade = reshape(fade,1,1,[]);
            ir = ir(:, :, 1:params.IRmax);
            HRTFs(iHRTF).SOFA.Data.IR = ir .* fade;
        end
        
        if ~isKey(HRTFGroupData, HRTFs(iHRTF).HRTFGroup)
            HRTFGroupData(HRTFs(iHRTF).HRTFGroup) = [];
        end
    
        HRTFGroupData(HRTFs(iHRTF).HRTFGroup) = [...
            HRTFGroupData(HRTFs(iHRTF).HRTFGroup) iHRTF];   
        
        fprintf('[%s][%s] azimuth: [%d, %d]; elevation: [%d, %d]; distance: %d\n', ...
            HRTFs(iHRTF).HRTFGroup, ...
            HRTFs(iHRTF).Name, ...
            min(HRTFs(iHRTF).Position(:, 1)), ...
            max(HRTFs(iHRTF).Position(:, 1)), ...
            min(HRTFs(iHRTF).Position(:, 2)), ...
            max(HRTFs(iHRTF).Position(:, 2)), ...
            HRTFs(iHRTF).Distance);
        
        if HRTFs(iHRTF).SOFA.Data.SamplingRate ~= params.RecordingsExpectedFs
            error('[%s][%s] Resampling from %d Hz to %d Hz', ...
                HRTF.HRTFGroup, HRTF.Name, ...
                HRTF.SOFA.Data.SamplingRate, ...
                params.RecordingsExpectedFs);
        end
    end
end


% Try load resampled HRTF routine
function [loadStatus, HRTF] = tryLoadResampledHRTF(id, HRTF, params)
    resampledSOFAdir = fullfile(params.HRTFBaseDir, ...
        ['_resampled_' num2str(params.RecordingsExpectedFs)], ...
        HRTF.HRTFGroup);
    resampledSOFAfilename = ['_resampled_' ...
        num2str(params.RecordingsExpectedFs) '_' HRTF.Name];
    fullSOFAfilename = fullfile(resampledSOFAdir, resampledSOFAfilename);
    
    if ~exist(fullSOFAfilename, 'file')
        loadStatus = false;
    else
        loadStatus = true;
        HRTF = loadHRTF(id, fullSOFAfilename, params);
    end
end


% Load HRTF routine
function HRTF = loadHRTF(id, filename, params)
    listing = dir(filename);
    fullFilename = fullfile(listing.folder, listing.name);
    filenameParts = split(listing.folder, filesep);
    SOFA = SOFAload(fullFilename);
    APV = SOFAcalculateAPV(SOFA);
    
    HRTF.Id = id;
    HRTF.Name = listing.name;
    HRTF.Folder = listing.folder;
    HRTF.HRTFGroup = filenameParts{end};
    HRTF.SOFA = SOFA;
    HRTF.Position = APV(:, 1:2);
    HRTF.Distance = unique(HRTF.SOFA.SourcePosition(:, 3));

    if any(strcmp(HRTF.HRTFGroup, params.InverseAzimuthHRTFGroups))
        HRTF.Position = HRTF.Position * [-1 0; 0 1];
    end
    
    if mod(HRTF.SOFA.API.N, 2) ~= 0
        tmpIR = HRTF.SOFA.Data.IR(:, :, 1:end-1); % Remove last sample
        HRTF.SOFA.Data.IR = tmpIR;
        HRTF.SOFA.API.N = size(tmpIR, 3);
    end    
end


% Resample and save routine
function HRTF = resampleAndSave(HRTF, params)    
    fprintf('[%s][%s] Resampling from %d Hz to %d Hz\n', ...
        HRTF.HRTFGroup, HRTF.Name, ...
        HRTF.SOFA.Data.SamplingRate, ...
        params.RecordingsExpectedFs);
    
    HRTF.SOFA = SOFAresample(HRTF.SOFA, params.RecordingsExpectedFs);
    
    resampledSOFAdir = fullfile(params.HRTFBaseDir, ...
        ['_resampled_' num2str(params.RecordingsExpectedFs)], ...
        HRTF.HRTFGroup);
    resampledSOFAfilename = ['_resampled_' ...
        num2str(params.RecordingsExpectedFs) '_' HRTF.Name];
    
    if ~exist(resampledSOFAdir, 'dir')
        mkdir(resampledSOFAdir);
    end
    
    fullSOFAfilename = fullfile(resampledSOFAdir, resampledSOFAfilename);
    HRTF.SOFA = SOFAsave(fullSOFAfilename, HRTF.SOFA, 0);
end


% Resample SOFA routine
function Obj = SOFAresample(Obj, targetFs)    
    currentFs = Obj.Data.SamplingRate;
    
    if currentFs == targetFs
        return
    end
    
    % Based on HRTFsamplingRateConverter10.m (S. Zieliński)
    M = size(Obj.Data.IR,1); % Number of measurements
    N = size(Obj.Data.IR,3); % Length of measurements
    IR = Obj.Data.IR;
    IR2 = zeros(M, 2, round(ceil(targetFs / currentFs * N*3) / 3));
    
    for ii = 1:M
        ir = squeeze(IR(ii, :, :))';
        irx3 = [ir; ir; ir];
        irx3 = resample(irx3, targetFs, currentFs);
        N2 = round(length(irx3)/3);
        ir2 = irx3(N2+1:2*N2, :);
        IR2(ii, :, :) = ir2';
    end
    
    Obj.Data.IR = IR2;
    Obj.Data.SamplingRate = targetFs;
    Obj=SOFAupdateDimensions(Obj);
end


% load tracks routine
function [tracks,trackNames] = loadAudioTracks(audioFilename, params)
    songName = fullfile(audioFilename.folder, audioFilename.name);
    trackFilenames = dir(fullfile(songName, '*.wav'));
%     disp(fullfile(trackFilenames(1).folder, trackFilenames(1).name))
    audioInfo = audioinfo(fullfile(trackFilenames(1).folder, ...
        trackFilenames(1).name));
    totalSamples = audioInfo.TotalSamples ...
        - params.RecordingLoadRange(1) * params.RecordingsExpectedFs;
    tracks = zeros(totalSamples, length(trackFilenames));
    
    for iTrackFilename = 1:length(trackFilenames)
        trackPath = fullfile(trackFilenames(iTrackFilename).folder, ...
            trackFilenames(iTrackFilename).name);
        [track,Fs] = audioread(trackPath, ...
            params.RecordingLoadRange * params.RecordingsExpectedFs + [1 0]);
        
        if Fs ~= params.RecordingsExpectedFs
            error('Track frequency is not expected frequency');
        end
        
        tracks(:, iTrackFilename) = track;
    end

    trackNames = {trackFilenames.name};
end

function tracks = normalizeAudioTracks(tracks, params)
    for iTrack = 1:size(tracks,2)
        tracks(:, iTrack) = normalizeLoudness(tracks(:, iTrack), ...
             params.RecordingsExpectedFs, ...
             params.TargetTrackLoudness);
    end
end

% Spatialize all audio trakcs routine
% spatResults shape (HRTF, dir, sample, ch)
function [outSpatResults,outSpatMetaResults] = spatializeSong(HRTFs, tracks, trackNames, audioName, params)    
    sz = [length(params.Elevations), ...
        length(HRTFs)];
    dur = params.RecordingSpatRange(2) * params.RecordingsExpectedFs;
    outSpatResults = zeros([sz ...
        dur params.NChannels]);
    outSpatMetaResults = cell(sz);
    trackNamesParts = split(trackNames, '.wav');
    trackNames = trackNamesParts(:, :, 1);
    
    for comb = allcomb(...
            1:length(params.Elevations), ...
            1:length(HRTFs))'
        cComb = num2cell(comb);
        [iElevation,iHRTF] = cComb{:};
        elevation = params.Elevations(iElevation);

        metaResults = getSceneMetaresult(HRTFs(iHRTF), audioName, ...
            trackNames, elevation, params);          
        
        spatResults = spatializeAudioTracks(...
            tracks, HRTFs(iHRTF), metaResults, params);

        outSpatResults(iElevation, iHRTF,:,:,:) = spatResults;
        outSpatMetaResults(iElevation, iHRTF) = {reshape(metaResults, 1, 1, [])};
        
        fprintf('Progress [el %d/%d] [HRTF %d/%d]\n', ...
            iElevation, ...
            length(params.Elevations), ...
            iHRTF, ...
            length(HRTFs));        
    end
    
    outSpatMetaResults = cell2mat(outSpatMetaResults);
end

function metaResults = getSceneMetaresult(HRTF, audioName, trackNames, elevation, params)
    width = rand * params.MaxWidth;    
    azimuthOffset = (params.AzimuthOffsetRange(2) ...
        - params.AzimuthOffsetRange(1)) .* rand ...
        + params.AzimuthOffsetRange(1);
    azimuthOffset = wrapTo180(azimuthOffset);

    randTrackAngles = rand(length(trackNames), 1);
    randTrackAngles = rescale(randTrackAngles, -width, width);
    randTrackAngles = randTrackAngles + azimuthOffset;
    randTrackAngles = wrapTo180(randTrackAngles);
    randTrackAngles(:, 2) = elevation;

    metaResults.AudioName = audioName;
    metaResults.TrackNames = trackNames;
    metaResults.HRTFId = HRTF.Id;
    metaResults.RandTrackAngles = randTrackAngles;
    metaResults.Elevation = elevation;
    metaResults.SceneWidth = width;
    metaResults.AzimuthEnsembleOffset = azimuthOffset;
end


% Spatialize audio routine
function spatResults = spatializeAudioTracks(...
    tracks, HRTF, metaResults, params)   

    spatResults = [];
    
    for iMetaresults = 1:length(metaResults)
        metaResult = metaResults(iMetaresults);
        interpHRTF = interpolateHRTF( ...
            HRTF.SOFA.Data.IR, ...
            HRTF.Position, ...
            metaResult.RandTrackAngles);
        spatResult = zeros(size(tracks,1) + size(interpHRTF,3) - 1, 2);        
        
        for iTrack = 1:size(interpHRTF, 1)
            track = tracks(:, iTrack); 
            spatTrack = [
                conv(squeeze(interpHRTF(iTrack, 1, :)), track) ...
                conv(squeeze(interpHRTF(iTrack, 2, :)), track)];
            spatResult = spatResult + spatTrack;
        end
        
        spatResult = trimAndFadeSignal(spatResult, params);

        if isempty(spatResults)
            spatResults = zeros([length(metaResults) size(spatResult)]);
        end
        spatResults(iMetaresults, :, :) = spatResult;
    end
end


% Trim and fade signal routine
function y = trimAndFadeSignal(x, params)    
    range = params.RecordingSpatRange * params.RecordingsExpectedFs - [0 1];
    y = x(range(1):sum(range), :);
    
    env = envGen(params.RecordingFadeTime(1), ...
        params.RecordingSpatRange(2), ...
        params.RecordingFadeTime(2), ...
        params.RecordingsExpectedFs, 2, 'sinsq')';
    y = y .* env;
end


% Postprocess spatialization results routine
function spatResults = posprocessSpatResults(spatResults, params)    
    % Peak normalization and scaling
    peakLevel = max(abs(spatResults), [], [3 4 5]);
    spatResults = params.RecordingLevelScale * spatResults ./ peakLevel;
    
    % DC equalization
    spatResults = spatResults - mean(spatResults, 3);
end


% Save spatialization results routine
% spatResults shape (HRTF, dir, sample, ch)
function spatResults = saveSpatResults(spatResults, spatMetaresults, ...
    HRTFs, audioFilename, iRep, params)
    
    if ~exist(params.SpatOutputDir, 'dir')
        mkdir(params.SpatOutputDir);
    end
    
    for comb = allcomb(...
            1:length(params.Elevations), ...
            1:length(HRTFs))'

        cComb = num2cell(comb);
        [iElevation,iHRTF] = cComb{:};
        elevation = params.Elevations(iElevation);        
        spatResult = squeeze(spatResults(iElevation, iHRTF, :, :));
        spatMeta = spatMetaresults(iElevation, iHRTF);
        acutalWidth = spatMeta.SceneWidth;
        azimuthOffset = spatMeta.AzimuthEnsembleOffset;
        
        spatFilename = getOutputFilename(...
            audioFilename, HRTFs(iHRTF), acutalWidth, ...
            elevation, azimuthOffset, iRep);
        
        spatParentDir = fullfile(params.SpatOutputDir);
        fullSpatFilename = fullfile(spatParentDir, spatFilename + '.wav');
        
        if ~exist(spatParentDir, 'dir')
            mkdir(spatParentDir);
        end       
        
        audiowrite(fullSpatFilename, spatResult, ...
            params.RecordingsExpectedFs, ...
            'BitsPerSample', 32);
    end
end

% Get output filename routine
function [filename,parentDir] = getOutputFilename( ...
    audioFilename, HRTF, ensembleActualWidth, elevation, azimuthOffset, ...
    iRep) 

    filename = sprintf("%s_hrtf%d_%s_width%.2f_el%.2f_azoffset%.2f_iteration%d", ...
        audioFilename.name, ...
        HRTF.Id, ...
        HRTF.HRTFGroup, ...
        ensembleActualWidth, ...
        elevation, ...
        azimuthOffset, ...
        iRep);
end


% Drawning routine
% spatMetaresults shape (width, HRTF, dir, audio)


function saveHRTFMetadata(HRTFs, params)
    T = struct2table(HRTFs);
    T = T(:, {'Id', 'Name', 'HRTFGroup', 'Distance'});
    filename = fullfile(params.FinalResultsOutputDir, 'HRTFs.csv');
    writetable(T, filename);
end

function signal = normalizeLoudness(signal, fs, targetLoudness)
    loudness = integratedLoudness(signal, fs);
    while abs(loudness - targetLoudness) > 0.001
        gain = 10^((targetLoudness - loudness)/20);
        signal = signal .* gain;
        loudness = integratedLoudness(signal, fs);
    end
end
