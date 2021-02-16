% %% Extract individual ERPs and plot over time
clear all; close all;

addpath(genpath('/Users/egrow1/Desktop/NEW_Redo_OriginalChans_SHAFTO/decode_libSVM'))
% Add the libSVM decoding library to the path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT: Settings and filenames for extracting ROI data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filepath = ['/Users/egrow1/Desktop/ALL_SW_noF_0_to_500ms_EBR/']; %set file path
%filepathResp = ['/Users/egrow1/Desktop/ALL_new_faceResp_SourceWaves/'];
hemisphere = 'Left'; % which hemisphere? Either 'Left' or 'Right', 'Both' (i.e all) or 'CombHemi' (combined each ROI bilaterally)

% Assign variables depending on which hemisphere being examined
if strcmp(hemisphere, 'Left')
    %Name the left/right coordinates to examine and labels and savename
    useCoords = {[-14,-69, 5]; [-32, -52, -6]; [-37, 41, 1]; [-10, 3, 66]}%; [-39, -84, 21]}; % Left hemisphere
    coordLabels = {'L_Occ','L_FFA', 'L_Inf_PFC', 'L_PFC'}%,'LITG','LPar','LOcc'};
     savename = ['AllEBResp_DiffVox_DecodeResults_SourcERP_Multifeature_' num2str(['L_Occ','L_FFA','L_InfPFC','L_SPPFC']) '.mat'];
elseif strcmp(hemisphere, 'Right')
    useCoords = {[33, -92, 10]; [30, -56, -6]; [8, 60, -22]; [17, -5, 61]}; % Left hemisphere
    coordLabels = {'R_Occ','R_FFA', 'R_Inf_PFC', 'R_PFC'}%,'LITG','LPar','LOcc'};
    savename = ['AllEBResp_DiffVox_DecodeResults_SourcERP_Multifeature_' num2str(['R_Occ','R_FFA','R_InfPFC','R_SPPFC']) '.mat'];
elseif strcmp(hemisphere, 'Both')
    useCoords = {[-14,-69, 5]; [-32, -52, -6]; [-22, 44, -13]; [-8, -6, 60]; ...
        [30, -96, 7]; [30, -56, -6]; [23, 45, -13]; [18, -6, 68]}; % Left hemisphere
    coordLabels = {'L_Occ','L_FFA', 'L_Inf_PFC', 'L_PFC',...
        'R_Occ','R_FFA', 'R_Inf_PFC', 'R_PFC'}%,'LITG','LPar','LOcc'};
    savename = ['AllEBResp_DiffVox_DecodeResults_SourcERP_Multifeature_AllROIs.mat'];
elseif strcmp(hemisphere, 'CombHemi')
    useCoords = {[-14, -69, 5],[30, -96, 7]; [-32, -52, -6],[30, -56, -6];...
        [-22, 44, -13], [23, 45, -13]; [-8, -6, 60], [18, -6, 68]}; % Left hemisphere
    coordLabels = {'L_Occ','L_FFA', 'L_Inf_PFC', 'L_PFC',...
        'R_Occ','R_FFA', 'R_Inf_PFC', 'R_PFC'};
    savename = ['AllEBResp_DiffVox_DecodeResults_SourcERP_Multifeature_AllBilateralROIsComb.mat'];
end

%Assign other parameters
zscore = 1; %1 = on, 0 = off
timeWindow = [0 600]; %examine this time window
timeIdx = linspace(0,600,301); % find index for time points in data
useWind = [find(timeIdx==timeWindow(1)):find(timeIdx==timeWindow(2))]; %assign idxs by time window
visualise = 0; %to show data visualisation =1, =0 for not
count = 1; %set up counter to plot (if visualise is turned on)

for gg = 1:length(useCoords)
    
    coordList = useCoords{gg}
    
    %% Do it for face trials
    for voxel = 1:size(coordList,1)
        thisVoxel = coordList;
        filename = ['Source_Waveform_at_' num2str(thisVoxel) '_for_ALL_trials_fmSae_faces_aEBfspmeeg_7MJN_0to600ms_trial_205.mat'];
        load ([filepath filename])
        GM_faces(voxel,:,:) = waveform;%(:,useWind); % time x trial x feature
    end
    
    %% Do it for RANDOM trials
    for voxel = 1:size(coordList,1)
        thisVoxel = coordList;
        filename = ['Source_Waveform_at_' num2str(thisVoxel) '_for_ALL_trials_fmSae_random_aEBfspmeeg_7MJN_0to600ms_trial_354.mat'];
        load ([filepath filename])
        GM_random(voxel,:,:) = waveform;%(:,useWind); % time x trial x feature
    end
    
    
    %% Gather all data within time window AND VISUALISE (if turned on)
    % GATHER DATA FOR DECODING
    GM_faces = squeeze(GM_faces(:,:,useWind)); %restrict to specified time window (if applicable)
    GM_random = squeeze(GM_random(:,:,useWind)); %restrict to specified time window (if applicable)
    
    data_faces = GM_faces'; %form features x trials
    data_random = GM_random'; %form features x trials
    
    if visualise == 1
        %% VISUALISE DATA ACROSS TIME AND TRIALS
        %VISUALISE within ENTIRE TIME WINDOW
        figure(1)
        subplot(3,4,count)
        set(gcf, 'Position', [0, 0, 1600, 1000])
        imagesc(data_faces',[-0.2 0.2]); hold on
        set(gca, 'FontSize', 20); colorbar;
        set(gca,'XTick',[1, 126, 251], 'XTickLabel', {'0', '250','500'});
        xlabel('Time (ms)'); ylabel('Voxel');
        thisTitle = [num2str(coordLabels{gg}) ' FACE TRIALS'];        
        title(thisTitle, 'FontSize',15); count = count + 1;
        subplot(3,4,count)
        imagesc(data_random',[-0.2 0.2]); hold on
        set(gca, 'FontSize', 20); colorbar;
        set(gca,'XTick',[1, 126, 251], 'XTickLabel', {'0', '250','500'});
        xlabel('Time (ms)'); ylabel('Voxel');
        thisTitle = [num2str(coordLabels{gg}) ' RANDOM TRIALS'];  
        title(thisTitle, 'FontSize',15); count = count + 1;
    end
    
    %% Decode usng libSVM
    nReps = 10;
    costRange = -5;
    
    [decodability, cv_acc, w] = decode_multiFeature_libSVM(data_faces, data_random, nReps, costRange, zscore)
  
    %Save all decoding data
    allDecode(gg,:) = decodability;
    allCV(gg,:,:) = cv_acc;
    allW{gg,:,:} = w;
    
    clearvars -except allW possCoords gg useCoords data_faces data_random filepath filePathCoords listCoords allDecode allCV ...
                filepathResp zscore useWind visualise coordLabels savename count saveFigDecode saveFigWaves
    
end

%% Plot these results 
%Find mean and std of decodability
decodability_mean = mean(allDecode')
decodability_std = std(allDecode')
decodability_sem = decodability_std/sqrt(10)

figure(2)
bar(decodability_mean); hold on;
ylim([45 100]); xlim([0 length(decodability_mean)+1])
errorbar(1:length(decodability_mean),decodability_mean, decodability_std,decodability_std)
set(gca, 'FontSize', 20); xlabel('Voxel ("TOP" single voxel)','FontSize',20); ylabel('Decoding Performance (%) (all T points, all trials)','FontSize',20);
set(gcf, 'Position', [0, 0, 1600, 1000])
plot(0:length(decodability_mean)+1, ones(1,length(decodability_mean)+2)*50, '-r', 'LineWidth', 2)
set(gca, 'XTick',[1:length(decodability_mean)], 'XTickLabel', coordLabels, 'FontSize',20);
title('Classification Accuracy - All Source Waveforms - Multifeature Decoding (all Time points, all Trials)', 'FontSize',20)

%Save these results
save(savename, 'decodability_mean', 'decodability_std', 'coordLabels', 'allCV','allDecode','allW', 'coordLabels')

%Plot using biolin plate
violinplot(allDecode',coordLabels)
  
