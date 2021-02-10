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

addpath(genpath('/Users/egrow1/Desktop/NEW_Redo_OriginalChans_SHAFTO/decode_libSVM'))
%Add the libSVM decoding library to the path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT: Settings and filenames for extracting ROI data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filepath = ['/Users/egrow1/Desktop/ALL_SW_8196_0_to_600ms_EBRem/']; %set file path
filepathList = ['/Users/egrow1/Desktop/NEW_Redo_OriginalChans_SHAFTO/SourceResults_MASSIVE/'];
filePathCoords = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/2_NonPara_Granger/'];
countA = 1; countB = 1; countC = 1;
maxFreq = 55; %maximum frequency to examine

% SELECT THE TWO NODES TO DECODE
possMainCoords = {'L_PFC','R_PFC','L_Inf_PFC','R_Inf_PFC'}%,'LeftPar','LeftOcc', ...
possSuppCoords = {'L_FFA','R_FFA','L_Occ','R_Occ'}; %'main' PFC coords

% % % % SELECT THE TWO NODES TO DECODE
% possMainCoords = {'R_Inf_PFC'}%,'LeftPar','LeftOcc', ...
% possSuppCoords = {'R_Occ'}; %'main' PFC coords

count = 1; cA = 1; cB = 1; cA_PFC = 1; cB_PFC = 1;
cC = 1; cD = 1; cC_PFC = 1; cD_PFC = 1;

for mainIdx = 1:length(possMainCoords)
    
    useMainCoord = possMainCoords{mainIdx}
    
    for suppIdx = 1:length(possSuppCoords)
        
        useSupCoord = possSuppCoords{suppIdx}

        % ASSIGN SAVE AND LOAD NAMES
         loadNameFaces = ['all_noF_faces_CohPower_GC_' num2str(useMainCoord) ...
        '_and_' num2str(useSupCoord) '_0to500ms_5Tapsmo_TopVoxel_FINAL.mat']; %savename
        loadNameRandom = ['all_noF_random_CohPower_GC_' num2str(useMainCoord) ...
        '_and_' num2str(useSupCoord) '_0to500ms_5Tapsmo_TopVoxel_FINAL.mat']; %savename

      savename = ['results_NoF_Decode_0to500ms_CohPowerGC_' num2str(useMainCoord) ...
                   '_and_' num2str(useSupCoord)  '_' num2str(maxFreq) 'Hz_5Tapsmo_TopVoxel_FINAL.mat']; %savename
        
        % ASSIGN DECODING PARAMETERS (libSVM)
        nReps = 10; % number of repetitions of decoding steps
        costRange = [-5]; %cost range = 1 value (%-20:5:15)
        zscore = 1; % zscore decoding data (1 = 'yes', 0 = 'no')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ANALYSIS: COHERENCE  - FACES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cohAll_faces = [];  cohAll_random = [];
        %% Load Face Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load([filePathCoords loadNameFaces]); %load face data
        
        %First look at coherence per trial
        freq_roi =  0 <= allGranger(1).freq & allGranger(1).freq < maxFreq
        freqWindow = allGranger(1).freq(freq_roi)
        
        for trialNo = 1:size(allCoherenceInfo,1)
            for thisFreq = 1:length(freqWindow)
                cohData_faces = allCoherenceInfo(trialNo).cohspctrm(:,:,thisFreq); %extract data
                cohAll_faces(trialNo,thisFreq,:) = cohData_faces(:)'; % coherence (form trials x freq x cohVals)
            end
        end
        
        cohAll_faces = cohAll_faces(:,:,2); % extract only value that matters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ANALYSIS: GC - FACES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        for trialNo = 1:size(allGranger,1)
            for thisFreq = 1:length(freqWindow)
                GCData_faces = allGranger(trialNo).grangerspctrm(:,:,thisFreq); % GC data
                GCAll_faces(trialNo,thisFreq,:) = GCData_faces(:)'; % Granger (form trials x freq x GCVals)
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ANALYSIS: POWER - FACES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        if mainIdx == 1 && suppIdx == 1 || mainIdx == 1 && suppIdx == 2 || mainIdx == 1 && suppIdx == 3 || ...
                mainIdx == 1 && suppIdx == 4 || mainIdx == 2 && suppIdx == 1 || mainIdx == 3 && suppIdx == 1 || ...
                mainIdx == 4 && suppIdx == 1
            
            for ss = 1:size(allPowerInfo,1)
                pwrData = allPowerInfo(ss).powspctrm(:,1:length(freqWindow))
                pwrAll_faces(ss,:,:) = pwrData'; % power (form trials x freq x pwrPerCond)
            end
            
            if mainIdx == 1 && suppIdx == 1 || mainIdx == 2 && suppIdx == 1 || ...
                    mainIdx == 3 && suppIdx == 1 || mainIdx == 4 && suppIdx == 1
                allROI_PowerPFC_faces(cA_PFC,:,:) = pwrAll_faces(:,:,2)';
                cA_PFC = cA_PFC + 1;
            end
            
            if mainIdx == 1 && suppIdx == 1 || mainIdx == 1 && suppIdx == 2 || ...
                    mainIdx == 1 && suppIdx == 3 || mainIdx == 1 && suppIdx == 4
                allROI_PowerROI_faces(cA,:,:) = pwrAll_faces(:,:,1)'; clear pwrAll_faces
                cA = cA + 1;
            end
        end
        
        allROI_coherence_faces(mainIdx,suppIdx,:,:) = cohAll_faces'; clear cohAll_faces
        allROI_Granger_toPFC_faces(mainIdx,suppIdx,:,:) = GCAll_faces(:,:,2)';
        allROI_Granger_fromPFC_faces(mainIdx,suppIdx,:,:) = GCAll_faces(:,:,3)'; clear GCAll_faces

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ANALYSIS: COHERENCE  - RANDOM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load([filePathCoords loadNameRandom]); %load random data
        
        
        for trialNo = 1:size(allCoherenceInfo,1)
            for thisFreq = 1:length(freqWindow)
                cohData_random = allCoherenceInfo(trialNo).cohspctrm(:,:,thisFreq)
                cohAll_random(trialNo,thisFreq,:) = cohData_random(:)'; % coherence (form trials x freq x cohVals)
            end
        end
        
        cohAll_random = cohAll_random(:,:,2); % extract only value that matters
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ANALYSIS: GC - FACES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
      
        for trialNo = 1:size(allGranger,1)
            for thisFreq = 1:length(freqWindow)
                GCData_random = allGranger(trialNo).grangerspctrm(:,:,thisFreq); % extract GC data
                GCAll_random(trialNo,thisFreq,:) = GCData_random(:)'; % Granger (form trials x freq x GCVals)
            end
        end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ANALYSIS: POWER - FACES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        if mainIdx == 1 && suppIdx == 1 || mainIdx == 1 && suppIdx == 2 || mainIdx == 1 && suppIdx == 3 || ...
                mainIdx == 1 && suppIdx == 4 || mainIdx == 2 && suppIdx == 1 || mainIdx == 3 && suppIdx == 1 || ...
                mainIdx == 4 && suppIdx == 1
            
            for ss = 1:size(allPowerInfo,1)
                pwrData_random = allPowerInfo(ss).powspctrm(:,1:length(freqWindow))
                pwrAll_random(ss,:,:) = pwrData_random'; % power (form trials x freq x pwrPerCond)
            end
            
            if mainIdx == 1 && suppIdx == 1 || mainIdx == 2 && suppIdx == 1 || ...
                    mainIdx == 3 && suppIdx == 1 || mainIdx == 4 && suppIdx == 1
                allROI_PowerPFC_random(cB_PFC,:,:) = pwrAll_random(:,:,2)';
                cB_PFC = cB_PFC + 1;
            end
            
            if mainIdx == 1 && suppIdx == 1 || mainIdx == 1 && suppIdx == 2 || ...
                    mainIdx == 1 && suppIdx == 3 || mainIdx == 1 && suppIdx == 4
                allROI_PowerROI_random(cB,:,:) = pwrAll_random(:,:,1)'; clear pwrAll_faces
                cB = cB + 1;
            end
        end
        
        allROI_coherence_random(mainIdx,suppIdx,:,:) = cohAll_random'; clear cohAll_random
        allROI_Granger_toPFC_random(mainIdx,suppIdx,:,:,:) = GCAll_random(:,:,2)';
        allROI_Granger_fromPFC_random(mainIdx,suppIdx,:,:,:) = GCAll_random(:,:,3)'; clear GCAll_random
        
        


%% FACES: SQUEEZE DATA INTO REQUIRED FORMAT (multiple ROIs)
% COHERENCE

% Gather COHERENCE data for both hemispheres
allCoh_faces = reshape(allROI_coherence_faces, [(size(allROI_coherence_faces,1))*(size(allROI_coherence_faces,2))*(size(allROI_coherence_faces,3)),size(allROI_coherence_faces,4)]);


% %POWER

% %Gather POWER for BOTH hemispheres
allPwr_atROI_faces = reshape(allROI_PowerROI_faces, [(size(allROI_PowerROI_faces,1))*(size(allROI_PowerROI_faces,2)),size(allROI_PowerROI_faces,3)]);
allLogPwr_atROI_faces = log10(allPwr_atROI_faces);

allPwr_atPFC_faces = reshape(allROI_PowerPFC_faces, [(size(allROI_PowerPFC_faces,1))*(size(allROI_PowerPFC_faces,2)),size(allROI_PowerPFC_faces,3)]);
allLogPwr_atPFC_faces = log10(allPwr_atPFC_faces);


%GRANGER

% %Gather GRANGER for BOTH hemispheres
allGC_toPFC_faces = reshape(allROI_Granger_toPFC_faces, [(size(allROI_Granger_toPFC_faces,1))*(size(allROI_Granger_toPFC_faces,2))*(size(allROI_Granger_toPFC_faces,3)),size(allROI_Granger_toPFC_faces,4)]);

allGC_fromPFC_faces = reshape(allROI_Granger_fromPFC_faces, [(size(allROI_Granger_fromPFC_faces,1))*(size(allROI_Granger_fromPFC_faces,2))*(size(allROI_Granger_fromPFC_faces,3)),size(allROI_Granger_fromPFC_faces,4)]);



%% RANDOM: SQUEEZE DATA INTO REQUIRED FORMAT
% Gather COHERENCE data for both hemispheres
allCoh_random = reshape(allROI_coherence_random, [(size(allROI_coherence_random,1))*(size(allROI_coherence_random,2))*(size(allROI_coherence_random,3)),size(allROI_coherence_random,4)]);


% %POWER

% %Gather POWER for BOTH hemispheres
allPwr_atROI_random = reshape(allROI_PowerROI_random, [(size(allROI_PowerROI_random,1))*(size(allROI_PowerROI_random,2)),size(allROI_PowerROI_random,3)]);
allLogPwr_atROI_random = log10(allPwr_atROI_random);

allPwr_atPFC_random = reshape(allROI_PowerPFC_random, [(size(allROI_PowerPFC_random,1))*(size(allROI_PowerPFC_random,2)),size(allROI_PowerPFC_random,3)]);
allLogPwr_atPFC_random = log10(allPwr_atPFC_random);


%GRANGER

% %Gather GRANGER for BOTH hemispheres
allGC_toPFC_random = reshape(allROI_Granger_toPFC_random, [(size(allROI_Granger_toPFC_random,1))*(size(allROI_Granger_toPFC_random,2))*(size(allROI_Granger_toPFC_random,3)),size(allROI_Granger_toPFC_random,4)]);

allGC_fromPFC_random = reshape(allROI_Granger_fromPFC_random, [(size(allROI_Granger_fromPFC_random,1))*(size(allROI_Granger_fromPFC_random,2))*(size(allROI_Granger_fromPFC_random,3)),size(allROI_Granger_fromPFC_random,4)]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS 1b: COHERENCE DECODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gather coherence data for decoding

% COHERENCE: FACE DATA FOR DECODING
cohDecode_faces = allCoh_faces; % = features (coh at freq) x trials
data.coh_faces = cohDecode_faces;

% COHERENCE: RANDOM DATA FOR DECODING
cohDecode_random = allCoh_random; % = features (coh at freq) x trials
data.coh_random = cohDecode_random;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COHERENCE: MULTI-FEATURE DECODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Decode COHERENCE usng libSVM
[decodability, cv_acc, w] = decode_multiFeature_libSVM(cohDecode_faces, cohDecode_random, nReps, costRange, zscore);

% Calculate mean decodability
decodability_mean = squeeze(mean(decodability)); % mean decoding of each channel after 10 repetitions
decodability_std = std(decodability); % std across each channel after repetitions
results.multiCoh = decodability;
results.multiCoh_mean = decodability_mean;
results.multiCoh_std = decodability_std;   %% Gather data (visualise if required)

cv_nu_mean = squeeze(mean(cv_acc,1))
cv_nu_std = std(cv_acc)

if length(costRange) > 1
    figure(11)
    errorbar(decodability_mean,decodability_std, 'LineWidth', 4); hold on;
    set(gca, 'FontSize', 30); xlabel('C (cost) parameter 2^x','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',20);
    set(gca, 'XTick',[1, 3, 5, 6], 'XTickLabel', {'-20','-10','0','5'}); xlim([1 length(decodability_mean)]);
    plot(1:length(decodability_mean)+1, ones(1,length(decodability_mean)+1)*50, '-r', 'LineWidth', 2); hold on;
    errorbar(cv_nu_mean,cv_nu_std, 'LineWidth', 4); hold on;
    title('COHERENCE: Varying Cost Parameter "C" for Combined Features', 'FontSize', 20); legend({'TEST', '', 'TRAIN'});
end


%% Visualise: Coherence multi-feature decoding results
figure(111)
bar(decodability_mean); hold on;
ylim([45 100]); xlim([0 length(decodability_mean)+1])
errorbar(1:length(decodability_mean),decodability_mean, decodability_std,decodability_std)
set(gca, 'FontSize', 30); xlabel('Coherence (multifeature)','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',20);
plot(0:length(decodability_mean)+1, ones(1,length(decodability_mean)+2)*50, '-r', 'LineWidth', 2)
title('COHERENCE Classification Accuracy - Face v Random Trials - Multifeature', 'FontSize', 20);

clear decodability; clear cv_acc;


%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS 2a: POWER RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS 2b: POWER DECODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gather power data for decoding

% POWER: FACE DATA FOR DECODING
pwr_NodeA_Decode_faces = allPwr_atROI_faces; %(:,:,1)'; % PWR at Node 1: features(power) x trials
pwr_NodeB_Decode_faces = allPwr_atPFC_faces;%(:,:,2)'; % PWR at Node 2: features(power) x trials
logPwr_NodeA_Decode_faces = allLogPwr_atROI_faces;%(:,:,1)'; % logPWR at Node 1: features(power) x trials
logPwr_NodeB_Decode_faces = allLogPwr_atPFC_faces;%(:,:,2)'; % logPWR at Node 2: features(power) x trials

%Save FACE data in structure
data.pwr_NodeA_faces = pwr_NodeA_Decode_faces;
data.pwr_NodeB_faces = pwr_NodeB_Decode_faces;
data.Logpwr_NodeA_faces = logPwr_NodeA_Decode_faces;
data.Logpwr_NodeB_faces = logPwr_NodeB_Decode_faces;

% POWER: RANDOM DATA FOR DECODING
pwr_NodeA_Decode_random = allPwr_atROI_random;%(:,:,1)'; %PWR at Node 1: features(power) x trials
pwr_NodeB_Decode_random = allPwr_atPFC_random;%(:,:,2)'; % PWR at Node 2: features(power) x trials
logPwr_NodeA_Decode_random = allLogPwr_atROI_random;%(:,:,1)'; % logPWR at Node 1: features(power) x trials
logPwr_NodeB_Decode_random = allLogPwr_atPFC_random;%(:,:,2)'; % logPWR at Node 2: features(power) x trials

%Save RANDOM data in structure
data.pwr_NodeA_random = pwr_NodeA_Decode_random;
data.pwr_NodeB_random = pwr_NodeB_Decode_random;
data.Logpwr_NodeA_random = logPwr_NodeA_Decode_random;
data.Logpwr_NodeB_random = logPwr_NodeB_Decode_random;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOG POWER: MULTI-FEATURE DECODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DECODE LOG POWER AT NODE 1
[decodability_NodeA, cv_acc] = decode_multiFeature_libSVM(logPwr_NodeA_Decode_faces, logPwr_NodeA_Decode_random, nReps, costRange, zscore);

% SAVE LOG POWER RESULTS AT NODE 1
decodability_mean_NodeA = squeeze(mean(decodability_NodeA)); % mean decoding of each channel after 10 repetitions
decodability_std_NodeA = std(decodability_NodeA); % std across each channel after repetitions

cv_nu_mean = squeeze(mean(cv_acc,1))
cv_nu_std = std(cv_acc)

if length(costRange) > 1
    figure(22)
    errorbar(decodability_mean_NodeA,decodability_std_NodeA, 'LineWidth', 4); hold on;
    set(gca, 'FontSize', 30); xlabel('C (cost) parameter 2^x','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',20);
    set(gca, 'XTick',[1, 3, 5], 'XTickLabel', {'-10','0','10'}); xlim([1 length(decodability_mean_NodeB_NodeA)]);
    plot(1:length(decodability_mean_NodeA)+1, ones(1,length(decodability_mean_NodeA)+1)*50, '-r', 'LineWidth', 2); hold on;
    errorbar(cv_nu_mean,cv_nu_std, 'LineWidth', 4); hold on;
    title('POWER Node A: Varying Cost Parameter "C" for Combined Features', 'FontSize', 20); legend({'TEST', '', 'TRAIN'});
end


%DECODE LOG POWER AT NODE 2
[decodability_NodeB, cv_acc] = decode_multiFeature_libSVM(logPwr_NodeB_Decode_faces, logPwr_NodeB_Decode_random, nReps, costRange, zscore);

% SAVE LOG POWER RESULTS AT NODE 2
decodability_mean_NodeB = squeeze(mean(decodability_NodeB)); % mean decoding of each channel after 10 repetitions
decodability_std_NodeB = std(decodability_NodeB); % std across each channel after repetitions

cv_nu_mean = squeeze(mean(cv_acc,1))
cv_nu_std = std(cv_acc)

if length(costRange) > 1
    figure(23)
    errorbar(decodability_mean_NodeB,decodability_std_NodeB, 'LineWidth', 4); hold on;
    set(gca, 'FontSize', 30); xlabel('C (cost) parameter 2^x','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',20);
    set(gca, 'XTick',[1, 3, 5], 'XTickLabel', {'-10','0','10'}); xlim([1 length(decodability_mean_NodeB_NodeA)]);
    plot(1:length(decodability_mean_NodeB)+1, ones(1,length(decodability_mean_NodeB)+1)*50, '-r', 'LineWidth', 2); hold on;
    errorbar(cv_nu_mean,cv_nu_std, 'LineWidth', 4); hold on;
    title('POWER Node B: Varying Cost Parameter "C" for Combined Features', 'FontSize', 20); legend({'TEST', '', 'TRAIN'});
end

% COMBINE POWER ACROSS ALL OF THE NODES
logPwr_NodeA_and_B_Decode_faces = [logPwr_NodeA_Decode_faces;logPwr_NodeB_Decode_faces]
logPwr_NodeA_and_B_Decode_random = [logPwr_NodeA_Decode_random;logPwr_NodeB_Decode_random]

%DECODE LOG POWER AT NODE 2
[decodability_NodeA_and_B, cv_acc] = decode_multiFeature_libSVM(logPwr_NodeA_and_B_Decode_faces, logPwr_NodeA_and_B_Decode_random, nReps, costRange, zscore);

% SAVE LOG POWER RESULTS AT NODE 2
decodability_mean_NodeA_and_B = squeeze(mean(decodability_NodeA_and_B)); % mean decoding of each channel after 10 repetitions
decodability_std_NodeA_and_B = std(decodability_NodeA_and_B); % std across each channel after repetitions


%% ALL RESULTS: SAVE THESE IN ONE STRUCTURE TOGETHER
decode_NodeA_NodeB = [decodability_mean_NodeA, decodability_mean_NodeB]
decode_NodeA_NodeB_std = [decodability_std_NodeA, decodability_std_NodeB]

%Save results in a structure
results.multiLogPwr_NodeA = decodability_NodeA;
results.multiLogPwr_NodeA_mean = decodability_mean_NodeA;
results.multiLogPwr_NodeA_std = decodability_std_NodeA;
results.multiLogPwr_NodeB = decodability_NodeB;
results.multiLogPwr_NodeB_mean = decodability_mean_NodeB;
results.multiLogPwr_NodeB_std = decodability_std_NodeB;

results.multiLogPwr_NodeA_and_B = decodability_NodeA_and_B;
results.multiLogPwr_NodeA_and_B_mean = decodability_mean_NodeA_and_B;
results.multiLogPwr_NodeA_and_B_std = decodability_std_NodeA_and_B;

%% Visualise: POWER - NODE 1 vs NODE 2
figure(2228)
bar(decode_NodeA_NodeB); hold on;
ylim([45 100]); xlim([0 length(decode_NodeA_NodeB)+1])
errorbar(1:length(decode_NodeA_NodeB),decode_NodeA_NodeB, decode_NodeA_NodeB_std,decode_NodeA_NodeB_std)
set(gca, 'FontSize', 30); xlabel('Power (multifeature)','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',20);
plot(0:length(decode_NodeA_NodeB_std)+1, ones(1,length(decode_NodeA_NodeB_std)+2)*50, '-r', 'LineWidth', 2); hold on;
set(gca, 'XTick', [1, 2], 'XTickLabel',{'FFA','PFC'});
title('LOG POWER Classification Accuracy - Face v Random Trials - Multifeature', 'FontSize', 20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS 3a: GRANGER CAUSALITY RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS 3b: GRANGER CAUSALITY DECODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gather GC data for decoding%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DATA GC FACES
GC_Decode_NodeB_to_NodeA_faces = allGC_fromPFC_faces;%(:,:,2)'; %features (GC at each freq) x trials
GC_Decode_NodeA_to_NodeB_faces = allGC_toPFC_faces;%(:,:,3)'; %features (GC at each freq) x trials

%DATA GC FACES: SAVE ALL THESE IN A STRUCTURE
data.GC_NodeB_to_NodeA_faces = GC_Decode_NodeB_to_NodeA_faces;
data.GC_NodeA_to_NodeB_faces = GC_Decode_NodeA_to_NodeB_faces;

% DATA GC RANDOM
GC_Decode_NodeB_to_NodeA_random = allGC_fromPFC_random;%(:,:,2)'; %features (GC at each freq) x trials
GC_Decode_NodeA_to_NodeB_random = allGC_toPFC_random; %(:,:,3)'; %features (GC at each freq) x trials

%DATA GC RANDOM: SAVE ALL THESE IN A STRUCTURE
data.GC_NodeB_to_NodeA_random = GC_Decode_NodeB_to_NodeA_random;
data.GC_NodeA_to_NodeB_random = GC_Decode_NodeA_to_NodeB_random;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRANGER CAUSALITY: MULTI FEATURE DECODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Decode GC at NODE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[decodability_NodeA_NodeB, cv_acc] = decode_multiFeature_libSVM(GC_Decode_NodeA_to_NodeB_faces, GC_Decode_NodeA_to_NodeB_random, nReps, costRange, zscore);

% SAVE THESE RESULTS
decodability_mean_NodeA_NodeB = squeeze(mean(decodability_NodeA_NodeB)); % mean decoding of each channel after 10 repetitions
decodability_std_NodeA_NodeB = std(decodability_NodeA_NodeB); % std across each channel after repetitions

cv_nu_mean = squeeze(mean(cv_acc,1))
cv_nu_std = std(cv_acc)

if length(costRange) > 1
    figure(33)
    errorbar(decodability_mean_NodeA_NodeB,decodability_std_NodeA_NodeB, 'LineWidth', 4); hold on;
    set(gca, 'FontSize', 30); xlabel('C (cost) parameter 2^x','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',20);
    set(gca, 'XTick',[1, 3, 5], 'XTickLabel', {'-10','0','10'}); xlim([1 length(decodability_mean_NodeB_NodeA)]);
    plot(1:length(decodability_mean_NodeA_NodeB)+1, ones(1,length(decodability_mean_NodeA_NodeB)+1)*50, '-r', 'LineWidth', 2); hold on;
    errorbar(cv_nu_mean,cv_nu_std, 'LineWidth', 4); hold on;
    title('GC Node A to Node B: Varying Cost Parameter "C" for Combined Features', 'FontSize', 20); legend({'TEST', '', 'TRAIN'});
end



%% Decode GC at NODE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[decodability_NodeB_NodeA, cv_acc] = decode_multiFeature_libSVM(GC_Decode_NodeB_to_NodeA_faces, GC_Decode_NodeB_to_NodeA_random, nReps, costRange, zscore);

%SAVE THESE RESULTS
decodability_mean_NodeB_NodeA = squeeze(mean(decodability_NodeB_NodeA)); % mean decoding of each channel after 10 repetitions
decodability_std_NodeB_NodeA = std(decodability_NodeB_NodeA); % std across each channel after repetitions

cv_nu_mean = squeeze(mean(cv_acc,1))
cv_nu_std = std(cv_acc)

if length(costRange) > 1
    figure(34)
    errorbar(decodability_mean_NodeB_NodeA,decodability_std_NodeB_NodeA, 'LineWidth', 4); hold on;
    set(gca, 'FontSize', 30); xlabel('C (cost) parameter 2^x','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',20);
    set(gca, 'XTick',[1, 3, 5], 'XTickLabel', {'-10','0','10'}); xlim([1 length(decodability_mean_NodeB_NodeA)]);
    plot(1:length(decodability_mean_NodeB_NodeA)+1, ones(1,length(decodability_mean_NodeB_NodeA)+1)*50, '-r', 'LineWidth', 2); hold on;
    errorbar(cv_nu_mean,cv_nu_std, 'LineWidth', 4); hold on;
    title('GC Node B to Node A: Varying Cost Parameter "C" for Combined Features', 'FontSize', 20); legend({'TEST', '', 'TRAIN'});
end

%SAVE ALL RESULTS TOGETHER
decode_GC_NodeA_NodeB = [decodability_mean_NodeA_NodeB, decodability_mean_NodeB_NodeA]
decode_GC_NodeA_NodeB_std = [decodability_std_NodeA_NodeB, decodability_std_NodeB_NodeA]

%SAVE THESE IN ONE STRUCTURE TOGETHER
results.multiGC_NodeA_to_NodeB = decodability_NodeA_NodeB;
results.multiGC_NodeA_to_NodeB_mean = decodability_mean_NodeA_NodeB;
results.multiGC_NodeA_to_NodeB_std = decodability_std_NodeA_NodeB;
results.multiGC_NodeB_to_NodeA = decodability_NodeB_NodeA;
results.multiGC_NodeB_to_NodeA_mean = decodability_mean_NodeB_NodeA;
results.multiGC_NodeB_to_NodeA_std = decodability_std_NodeB_NodeA;

%% Visualise: Multi-feature decoding Face v Random using GC TO/FROM NODE 1/2
figure(3333)
bar(decode_GC_NodeA_NodeB); hold on;
ylim([45 100]); xlim([0 length(decode_GC_NodeA_NodeB)+1])
errorbar(1:length(decode_GC_NodeA_NodeB),decode_GC_NodeA_NodeB, decode_GC_NodeA_NodeB_std,decode_GC_NodeA_NodeB_std)
set(gca, 'FontSize', 30); xlabel('GC (multifeature)','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',20);
plot(0:length(decode_GC_NodeA_NodeB_std)+1, ones(1,length(decode_GC_NodeA_NodeB_std)+2)*50, '-r', 'LineWidth', 2); hold on;
set(gca, 'XTick', [1, 2], 'XTickLabel',{'FFA to PFC','PFC to FFA'});
title('GC Classification Accuracy - Face v Random Trials - Multifeature', 'FontSize', 20);

%
%% SAVE DATA AND RESULTS
save(savename, 'results','data')
%close all;

    end
end
