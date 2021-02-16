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

% SELECT THE TWO NODES TO DECODE
possMainCoords = {'L_Inf_PFC','R_Inf_PFC','L_PFC','R_PFC'}; %'main' PFC coords
possSuppCoords = {'L_FFA','R_FFA','L_Occ','R_Occ'}%,'LeftPar','LeftOcc', ...

count = 1;

for mainIdx = 1:length(possMainCoords)
    
    useMainCoord = possMainCoords{mainIdx}
    
    for suppIdx = 1:length(possSuppCoords)
        
        useSupCoord = possSuppCoords{suppIdx}

        % ASSIGN SAVE AND LOAD NAMES
         loadNameFaces = ['allEB_DiffVoxels_faces_CohPower_GC_' num2str(useMainCoord) ...
        '_and_' num2str(useSupCoord) '_TopVoxel_FINAL.mat']; %savename
        loadNameRandom = ['allEB_DiffVoxels_random_CohPower_GC_' num2str(useMainCoord) ...
        '_and_' num2str(useSupCoord) '_TopVoxel_FINAL.mat']; %savename

    savename = ['results_REAL_v_PERMUTATION_Decode_CohPowerGC_CoordsALL_MultiROI_TopVoxel_FINAL.mat']; %savename
        
        % ASSIGN DECODING PARAMETERS (libSVM)
        nReps = 10; % number of repetitions of decoding steps
        costRange = [-5]; %cost range = 1 value (%-20:5:15)
        zscore = 1; % zscore decoding data (1 = 'yes', 0 = 'no')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ANALYSIS 1a: COHERENCE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cohAll_faces = [];  cohAll_random = [];
        %% Load Face Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load([filePathCoords loadNameFaces]); %load face data
        
        %First look at coherence per trial
        freq_roi =  0 <= allGranger(1).freq & allGranger(1).freq < 40
        freqWindow = allGranger(1).freq(freq_roi)
        
        for trialNo = 1:size(allCoherenceInfo,1)
            for thisFreq = 1:length(freqWindow)
                cohData_faces = allCoherenceInfo(trialNo).cohspctrm(:,:,thisFreq); %extract data
                cohAll_faces(trialNo,thisFreq,:) = cohData_faces(:)'; % coherence (form trials x freq x cohVals)
            end
        end
        
        cohAll_faces = cohAll_faces(:,:,2); % extract only value that matters
        
        for trialNo = 1:size(allGranger,1)
            for thisFreq = 1:length(freqWindow)
                GCData_faces = allGranger(trialNo).grangerspctrm(:,:,thisFreq); % GC data
                GCAll_faces(trialNo,thisFreq,:) = GCData_faces(:)'; % Granger (form trials x freq x GCVals)
            end
        end
        
        
        for ss = 1:size(allPowerInfo,1)
            pwrData = allPowerInfo(ss).powspctrm(:,1:length(freqWindow))
            pwrAll_faces(ss,:,:) = pwrData'; % power (form trials x freq x pwrPerCond)
        end
        
        allROI_coherence_faces(mainIdx,suppIdx,:,:) = cohAll_faces'; clear cohAll_faces
        allROI_Granger_toPFC_faces(mainIdx,suppIdx,:,:) = GCAll_faces(:,:,2)';
        allROI_Granger_fromPFC_faces(mainIdx,suppIdx,:,:) = GCAll_faces(:,:,3)'; clear GCAll_faces
        allROI_PowerPFC_faces(mainIdx,suppIdx,:,:) = pwrAll_faces(:,:,2)';
        allROI_PowerROI_faces(mainIdx,suppIdx,:,:) = pwrAll_faces(:,:,1)'; clear pwrAll_faces
        
        
        %% GATHER RESULTS FOR RANDOM DATA
        load([filePathCoords loadNameRandom]); %load random data
        
        
        for trialNo = 1:size(allCoherenceInfo,1)
            for thisFreq = 1:length(freqWindow)
                cohData_random = allCoherenceInfo(trialNo).cohspctrm(:,:,thisFreq)
                cohAll_random(trialNo,thisFreq,:) = cohData_random(:)'; % coherence (form trials x freq x cohVals)
            end
        end
        
        cohAll_random = cohAll_random(:,:,2); % extract only value that matters
        
        
        for trialNo = 1:size(allGranger,1)
            for thisFreq = 1:length(freqWindow)
                GCData_random = allGranger(trialNo).grangerspctrm(:,:,thisFreq); % extract GC data
                GCAll_random(trialNo,thisFreq,:) = GCData_random(:)'; % Granger (form trials x freq x GCVals)
            end
        end
        
        for ss = 1:size(allPowerInfo,1)
            pwrData_random = allPowerInfo(ss).powspctrm(:,1:length(freqWindow))
            pwrAll_random(ss,:,:) = pwrData_random'; % power (form trials x freq x pwrPerCond)
        end
        
        allROI_coherence_random(mainIdx,suppIdx,:,:) = cohAll_random'; clear cohAll_random
        allROI_Granger_toPFC_random(mainIdx,suppIdx,:,:,:) = GCAll_random(:,:,2)';
        allROI_Granger_fromPFC_random(mainIdx,suppIdx,:,:,:) = GCAll_random(:,:,3)'; clear GCAll_random
        allROI_PowerPFC_random(mainIdx,suppIdx,:,:) = pwrAll_random(:,:,2)';
        allROI_PowerROI_random(mainIdx,suppIdx,:,:) = pwrAll_random(:,:,1)'; clear pwrAll_random
        
        
        
    end
end

nPerms = 100;

for kk = 1:nPerms
%% PERMUTATION ANALYSIS: ORIGINAL DATA FACES
% COHERENCE
% Gather COHERENCE data for both hemispheres
allCoh_faces = reshape(allROI_coherence_faces, [(size(allROI_coherence_faces,1))*(size(allROI_coherence_faces,2))*(size(allROI_coherence_faces,3)),size(allROI_coherence_faces,4)]);

% %POWER
% %Gather POWER for BOTH hemispheres
allPwr_atROI_faces = reshape(allROI_PowerROI_faces, [(size(allROI_PowerROI_faces,1))*(size(allROI_PowerROI_faces,2))*(size(allROI_PowerROI_faces,3)),size(allROI_PowerROI_faces,4)]);
allLogPwr_atROI_faces = log10(allPwr_atROI_faces);

allPwr_atPFC_faces = reshape(allROI_PowerPFC_faces, [(size(allROI_PowerPFC_faces,1))*(size(allROI_PowerPFC_faces,2))*(size(allROI_PowerPFC_faces,3)),size(allROI_PowerPFC_faces,4)]);
allLogPwr_atPFC_faces = log10(allPwr_atPFC_faces);


%GRANGER
% %Gather GRANGER for BOTH hemispheres
allGC_toPFC_faces = reshape(allROI_Granger_toPFC_faces, [(size(allROI_Granger_toPFC_faces,1))*(size(allROI_Granger_toPFC_faces,2))*(size(allROI_Granger_toPFC_faces,3)),size(allROI_Granger_toPFC_faces,4)]);
allGC_fromPFC_faces = reshape(allROI_Granger_fromPFC_faces, [(size(allROI_Granger_fromPFC_faces,1))*(size(allROI_Granger_fromPFC_faces,2))*(size(allROI_Granger_fromPFC_faces,3)),size(allROI_Granger_fromPFC_faces,4)]);


%% PERMUTATION ANALYSIS: ORIGINAL DATA RANDOM
% Gather COHERENCE data for both hemispheres
allCoh_random = reshape(allROI_coherence_random, [(size(allROI_coherence_random,1))*(size(allROI_coherence_random,2))*(size(allROI_coherence_random,3)),size(allROI_coherence_random,4)]);

% %POWER
% %Gather POWER for BOTH hemispheres
allPwr_atROI_random = reshape(allROI_PowerROI_random, [(size(allROI_PowerROI_random,1))*(size(allROI_PowerROI_random,2))*(size(allROI_PowerROI_random,3)),size(allROI_PowerROI_random,4)]);
allLogPwr_atROI_random = log10(allPwr_atROI_random);

allPwr_atPFC_random = reshape(allROI_PowerPFC_random, [(size(allROI_PowerPFC_random,1))*(size(allROI_PowerPFC_random,2))*(size(allROI_PowerPFC_random,3)),size(allROI_PowerPFC_random,4)]);
allLogPwr_atPFC_random = log10(allPwr_atPFC_random);


%GRANGER

% %Gather GRANGER for BOTH hemispheres
allGC_toPFC_random = reshape(allROI_Granger_toPFC_random, [(size(allROI_Granger_toPFC_random,1))*(size(allROI_Granger_toPFC_random,2))*(size(allROI_Granger_toPFC_random,3)),size(allROI_Granger_toPFC_random,4)]);

allGC_fromPFC_random = reshape(allROI_Granger_fromPFC_random, [(size(allROI_Granger_fromPFC_random,1))*(size(allROI_Granger_fromPFC_random,2))*(size(allROI_Granger_fromPFC_random,3)),size(allROI_Granger_fromPFC_random,4)]);


%% PERMUTATION ANALYSIS: PERM DATA FACES and RANDOM

% PERMUTE THE COHERENCE DATA
allCoh = [allCoh_faces, allCoh_random];

randPerm = reshape(randperm(length(allCoh(:))), size(allCoh,1), size(allCoh,2));
permCoh = allCoh(randPerm);

permCoh_faces = permCoh(1:size(allCoh_faces,1), 1:size(allCoh_faces,2));
permCoh_random = permCoh(1:size(allCoh_random,1), size(allCoh_faces,2)+1: end);

%PERMUTE THE LOG POWER DATA AT ROI
allLowPwr_ROI = [allLogPwr_atROI_faces, allLogPwr_atROI_random];

randPermPwrROI = reshape(randperm(length(allLowPwr_ROI(:))), size(allLowPwr_ROI,1), size(allLowPwr_ROI,2));
permPwrROI = allLowPwr_ROI(randPermPwrROI);

permLogPwr_ROI_faces = permCoh(1:size(allLogPwr_atROI_faces,1), 1:size(allLogPwr_atROI_faces,2));
permLogPwr_ROI_random = permCoh(1:size(allLogPwr_atROI_random,1), size(allLogPwr_atROI_faces,2)+1: end);

%PERMUTE THE LOG POWER DATA AT PFC
allLowPwr_PFC = [allLogPwr_atPFC_faces, allLogPwr_atPFC_random];

randPermPwrPFC = reshape(randperm(length(allLowPwr_PFC(:))), size(allLowPwr_PFC,1), size(allLowPwr_PFC,2));
permPwrPFC = allLowPwr_PFC(randPermPwrPFC);

permLogPwr_PFC_faces = permCoh(1:size(allLogPwr_atPFC_faces,1), 1:size(allLogPwr_atPFC_faces,2));
permLogPwr_PFC_random = permCoh(1:size(allLogPwr_atPFC_random,1), size(allLogPwr_atPFC_faces,2)+1: end);

%PERMUTE THE GC ESTIMATE DATA TO PFC
allGC_toPFC = [allGC_toPFC_faces, allGC_toPFC_random];

randGCtoPFC = reshape(randperm(length(allGC_toPFC(:))), size(allGC_toPFC,1), size(allGC_toPFC,2));
permGCtoPFC = allGC_toPFC(randGCtoPFC);

permGCto_PFC_faces = permCoh(1:size(allGC_toPFC_faces,1), 1:size(allGC_toPFC_faces,2));
permGCto_PFC_random = permCoh(1:size(allGC_toPFC_random,1), size(allGC_toPFC_faces,2)+1: end);


%PERMUTE THE GC ESTIMATE DATA FROM PFC
allGC_fromPFC = [allGC_fromPFC_faces, allGC_fromPFC_random];

randGCfromPFC = reshape(randperm(length(allGC_fromPFC(:))), size(allGC_fromPFC,1), size(allGC_fromPFC,2));
permGCfromPFC = allGC_fromPFC(randGCfromPFC);

permGCfrom_PFC_faces = permCoh(1:size(allGC_fromPFC_faces,1), 1:size(allGC_fromPFC_faces,2));
permGCfrom_PFC_random = permCoh(1:size(allGC_fromPFC_random,1), size(allGC_fromPFC_faces,2)+1: end);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS 1b: COHERENCE DECODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COHERENCE: MULTI-FEATURE DECODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Decode COHERENCE usng libSVM
[decodability_coh] = decode_multiFeature_libSVM(allCoh_faces, allCoh_random, nReps, costRange, zscore);
mean_decode_coh(kk) = mean(decodability_coh)

[decodability_perm_coh] = decode_multiFeature_libSVM(permCoh_faces, permCoh_random, nReps, costRange, zscore);
mean_decode_perm_coh(kk) = mean(decodability_perm_coh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS 2b: LOG POWER DECODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DECODE LOG POWER AT NODE 1
[decodability_pwr_ROI] = decode_multiFeature_libSVM(allLogPwr_atROI_faces, allLogPwr_atROI_random, nReps, costRange, zscore);
mean_decode_pwr_ROI(kk) = mean(decodability_pwr_ROI)

[decodability_perm_pwr_ROI] = decode_multiFeature_libSVM(permLogPwr_ROI_faces, permLogPwr_ROI_random, nReps, costRange, zscore);
mean_decode_perm_pwr_ROI(kk) = mean(decodability_perm_pwr_ROI)

%DECODE LOG POWER AT NODE 2
[decodability_pwr_PFC] = decode_multiFeature_libSVM(allLogPwr_atPFC_faces, allLogPwr_atPFC_random, nReps, costRange, zscore);
mean_decode_pwr_PFC(kk) = mean(decodability_pwr_PFC)

[decodability_perm_pwr_PFC] = decode_multiFeature_libSVM(permLogPwr_PFC_faces, permLogPwr_PFC_random, nReps, costRange, zscore);
mean_decode_perm_pwr_PFC(kk) = mean(decodability_perm_pwr_PFC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS 3a: GRANGER CAUSALITY RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Decode GC at NODE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[decodability_GC_toPFC] = decode_multiFeature_libSVM(allGC_toPFC_faces, allGC_toPFC_random, nReps, costRange, zscore);
mean_decode_GC_toPFC(kk) = mean(decodability_GC_toPFC)

[decodability_perm_GC_toPFC] = decode_multiFeature_libSVM(permGCto_PFC_faces, permGCto_PFC_random, nReps, costRange, zscore);
mean_decode_perm_GC_toPFC(kk) = mean(decodability_perm_GC_toPFC)

%% Decode GC at NODE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[decodability_GC_fromPFC] = decode_multiFeature_libSVM(allGC_fromPFC_faces, allGC_fromPFC_random, nReps, costRange, zscore);
mean_decode_fromPFC(kk) = mean(decodability_GC_fromPFC)

[decodability_perm_GC_fromPFC] = decode_multiFeature_libSVM(permGCfrom_PFC_faces, permGCfrom_PFC_random, nReps, costRange, zscore);
mean_decode_perm_GC_fromPFC(kk) = mean(decodability_perm_GC_fromPFC)



end
%


%% TEST SIGNIFICANCE OF REAL versus PERMUTED RESULTS
% COHERENCE
[H_coh,P_coh,CI_coh,STATS_coh] = ttest2(mean_decode_coh, mean_decode_perm_coh)
final_P_coh = P_coh*nPerms

% LOG PWR PFC
[H_pwrPFC,P_pwrPFC,CI_pwrPFC,STATS_pwrPFC] = ttest2(mean_decode_pwr_PFC, mean_decode_perm_pwr_PFC)
final_P_pwrPFC = P_pwrPFC*nPerms

% LOG PWR ROIS
[H_pwrROI,P_pwrROI,CI_pwrROI,STATS_pwrROI] = ttest2(mean_decode_pwr_ROI, mean_decode_perm_pwr_ROI)
final_P_pwrROI = P_pwrROI*nPerms

% GC TO PFC
[H_GCtoPFC,P_GCtoPFC,CI_GCtoPFC,STATS_GCtoPFC] = ttest2(mean_decode_GC_toPFC, mean_decode_perm_GC_toPFC)
final_P_GCtoPFC = P_GCtoPFC*nPerms

% GC FROM PFC

[H_GCfromPFC,P_GCfromPFC,CI_GCfromPFC,STATS_GCfromPFC] = ttest2(mean_decode_fromPFC, mean_decode_perm_GC_fromPFC)
final_P_GCfromPFC = P_GCfromPFC*nPerms


%% SAVE DATA AND RESULTS
save(savename)
%close all;

