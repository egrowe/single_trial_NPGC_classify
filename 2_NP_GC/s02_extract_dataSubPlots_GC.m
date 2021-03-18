%% Visualise and decode NP coherence, power and Granger Causality results

% This script enable visualisation of the results from non-parametric
% Granger Causality analysis including coherence, power and GC
% coefficients.

% This script requires the libSVM toolbox is added to the path

% Each stage split into VISUALISATION raw data, DECODING data and
% visualisation of the decoding results

% Extract data for visualisation and decoding for either:
% (1): Coherence
% (2): Log Power
% (3): Granger causality

close all; clear all;

%Add the libSVM decoding library to the path
addpath(genpath('/Users/egrow1/Desktop/NEW_Redo_OriginalChans_SHAFTO/decode_libSVM'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT: Settings and filenames for extracting ROI data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filepath = ['/Users/egrow1/Desktop/GC_estimates_OSF/']; %set file path
filepathList = ['/Users/egrow1/Desktop/NEW_Redo_OriginalChans_SHAFTO/SourceResults_MASSIVE/'];
filePathCoords = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/1_EEG_PreProcessing/'];

% SELECT THE TWO NODES TO DECODE GC ESTIMATES
possSensoryCoords = {'L_Occ','R_Occ','L_FFA', 'R_FFA'} %'main' PFC coords
possPFCCoords = {'L_Inf_PFC','L_PFC','R_Inf_PFC','R_PFC'}; %,'LeftPar','LeftOcc', ...
maxFreq = 55; %maximum frequency to examine 

%Setup counters
countA = 1; countB = 1; countC = 1;

for mainIdx = 1:length(possSensoryCoords)
    useSensoryCoord = possSensoryCoords{mainIdx}; %select this sensory ROI location
    for suppIdx = 1:length(possPFCCoords)
        usePFCCoord = possPFCCoords{suppIdx}; %select this PFC ROI location
        
        % ASSIGN LOAD AND SAVE NAMES
        loadNameFaces = ['GCest_faces_' num2str(usePFCCoord) ...
            '_and_' num2str(useSensoryCoord) '_0to500ms_4Tapsmo_TopVoxel.mat']; %load faces
        loadNameRandom = ['GCest_random_' num2str(usePFCCoord) ...
            '_and_' num2str(useSensoryCoord) '_0to500ms_4Tapsmo_TopVoxel.mat']; %load random
        
        savename = ['resultsDecode_' num2str(usePFCCoord) ...
            '_and_' num2str(useSensoryCoord) '_0to500ms_4Tapsmo_TopVoxel.mat']; %savename
                
        % ASSIGN DECODING PARAMETERS (libSVM)
        nReps = 10; % number of repetitions of decoding steps
        costRange = [-5]; %cost range = 1 value (%-20:5:15) set to range if doing optimisation
        zscore = 1; % zscore decoding data (1 = 'yes', 0 = 'no')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ANALYSIS 1a: GRANGER CAUSALITY RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Load Face Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load([filepath loadNameFaces]); %load face data
        
        %First determine the frequency window to examine
        freq_roi =  0 <= allGranger(1).freq & allGranger(1).freq < maxFreq
        freqWindow = allGranger(1).freq(freq_roi)

        for trialNo = 1:size(allGranger,1)
            for thisFreq = 1:length(freqWindow)
                GCData_faces = allGranger(trialNo).grangerspctrm(:,:,thisFreq); % GC data
                GCAll_faces(trialNo,thisFreq,:) = GCData_faces(:)'; % Granger (form trials x freq x GCVals)
            end
        end
        
        %Gather mean and std of GC
        mGC_faces = squeeze(mean(GCAll_faces,1)); %MEAN GC over trials (at each frequency band)
        stdGC_faces = squeeze(std(GCAll_faces)); %STD GC over trials
        semGC_faces = stdGC_faces/(sqrt(length(GCAll_faces))); %SEM GC over trials

        GC_TO_PFC_faces = mGC_faces(:,3);
        GC_sem_TO_PFC_faces = semGC_faces(:,3);
        GC_FROM_PFC_faces = mGC_faces(:,2);
        GC_sem_FROM_PFC_faces = semGC_faces(:,2);
        
        figure(333)
        %thisPlot = usePlotIdxs(countA);
        subplot(6,4,countC)
        errorbar(GC_TO_PFC_faces,GC_sem_TO_PFC_faces, 'LineWidth', 3); hold on; ylim([0.2 0.75]);
        if suppIdx == 1 && mainIdx == 1
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {ceil(freqWindow(1)),ceil(freqWindow(10)),ceil(freqWindow(20))}); xlim([1 length(mCoh_faces)]);
        else set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]); xlim([1 length(mCoh_faces)]);
        end
        % title(['GC to ' num2str(useMainCoord) ' from:' num2str(useSupCoord)], 'FontSize', 15); %legend({'', 'TEST', '', 'TRAIN'});
        hold on;
        
        figure(444)
        subplot(6,4,countC)
        errorbar(GC_FROM_PFC_faces,GC_sem_FROM_PFC_faces, 'LineWidth', 3); hold on; ylim([0.2 0.75]);
        if suppIdx == 1 && mainIdx == 1
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {ceil(freqWindow(1)),ceil(freqWindow(10)),ceil(freqWindow(20))}); xlim([1 length(mCoh_faces)]);
        else set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]); xlim([1 length(mCoh_faces)]);
        end
        %  title(['GC from ' num2str(useMainCoord) ' to: ' num2str(useSupCoord)], 'FontSize', 15); %legend({'', 'TEST', '', 'TRAIN'});
        hold on;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Load Random Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load([filepath loadNameRandom]); %load random data
        
        for trialNo = 1:size(allGranger,1)
            for thisFreq = 1:length(freqWindow)
                GCData_random = allGranger(trialNo).grangerspctrm(:,:,thisFreq); % extract GC data
                GCAll_random(trialNo,thisFreq,:) = GCData_random(:)'; % Granger (form trials x freq x GCVals)
            end
        end
        
        %Gather mean and std of GC
        mGC_random = squeeze(mean(GCAll_random,1))
        stdGC_random = squeeze(std(GCAll_random))
        semGC_random = stdGC_random/(sqrt(length(GCAll_random)))

        %Gather GC data
        GC_TO_PFC_random = mGC_random(:,3);
        GC_sem_TO_PFC_random = semGC_random(:,3);
        GO_FROM_PFC_random = mGC_random(:,2);
        GC_sem_FROM_PFC_random = semGC_random(:,2);
        
        figure(333)
        subplot(6,4,countC)
        errorbar(GC_TO_PFC_random,GC_sem_TO_PFC_random, 'LineWidth', 3); hold on;
        
        figure(444)
        subplot(6,4,countC)
        errorbar(GO_FROM_PFC_random,GC_sem_FROM_PFC_random, 'LineWidth', 3); hold on;
        
        countC = countC + 1;
        
         end
end
