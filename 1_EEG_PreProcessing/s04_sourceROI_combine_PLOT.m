%% Determine which voxels have the highest t-stat value over time
%   (use this as the representative voxel for the ROI -- used for GC and
%   decodng)

% NOTE: Before executing this script, all 8196 voxel coordinates need to be
% assigned to an ROI. For this, use the batch coordinate input option of
% the AnatomyToolbox extension for SPM12. The output of this file will be
% the assigned cortical locations for each voxel.

% Next, manually sort these coordinates to ROIs. In this study, we selected
% Left and Right PFC (defined as posterior, superiorl, medial prefrontal
% without inferior or orbital prefrontal), FFA, ITG, parietal and
% occipital. Once manually sorted, the coordinates for each ROI need to be
% saved as separate .mat files that are loaded in this script to determine
% the 'top voxel' per ROI (defined as the voxel with the highest t-stat
% value when cmoparing face and random trials over the 0 to 500 ms time
% window -- top t-stat at any point within this time window).

clear all; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT: Settings and filenames for extracting ROI data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filepath = ['/Users/egrow1/Desktop/ALL_SW_noF_0_to_500ms_EBR/']; %set file path
filepathCoords = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/1_EEG_PreProcessing/'];
listCoords = ['Coordinates_FINAL_UPDATED_8196_voxels_byROI.mat'];
load([filepathCoords listCoords])

%Extract coordinates for these ROIs
%coordROI = {'R_PFC','R_Inf_PFC','R_FFA','R_Occ'}
coordROI = {'L_PFC','L_Inf_PFC','L_FFA','L_Occ'}


%Define other parameters
zscore = 1; %1 = on, 0 = offclea
timeWindow = [0 500]; %examine this time window
timeIdx = linspace(0,500,151); % find index for time points in data
useWind = [find(timeIdx==timeWindow(1)):find(timeIdx==timeWindow(2))]; %assign idxs by time window
visualise = 1; %to show data visualisation =1, =0 for not

for ff = 1:length(coordROI)
    
    useCoordName = coordROI{ff}
    load(['Coord_NoF_FINAL_' num2str(useCoordName) '_by_TstatOverTime_0_to_500ms.mat'], 'sortedCoords');
    thisCoord = sortedCoords(1,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Extract waveforms for face trials
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        thisVoxel = thisCoord;
        filename = ['NEW_Source_Waveform_0to500ms_at_' num2str(thisVoxel) '_for_ALL_trials_mSae_faces_aEBfspmeeg_7MJN_0to500ms_trial_205.mat'];
        load ([filepath filename])
        GM_faces(ff,:,:) = waveform; % time x trial x feature
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Extract waveforms for RANDOM trials
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        filename = ['NEW_Source_Waveform_0to500ms_at_' num2str(thisVoxel) '_for_ALL_trials_mSae_random_aEBfspmeeg_7MJN_0to500ms_trial_354.mat'];
        load ([filepath filename])
        GM_random(ff,:,:) = waveform; % time x trial x feature
    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gather all data within time window AND VISUALISE

%Set colours for the plots
faces_color = [0.3, 0.1, 0.3]; %purple
random_color = [0.1, 0.4, 0.4]; %green


for gg = 1:length(coordROI)
    
    faceData = squeeze(GM_faces(gg,:,:)); %extract face data for this coord
    randData =  squeeze(GM_random(gg,:,:)); %extract random data this coord
    
    %Plot using JackKnife plot
    subplot(1,length(coordROI),gg)
    plot_JackKnife(faceData, faces_color, '-', 0:2:500)
    hold on;
    plot_JackKnife(randData, random_color, '-', 0:2:500)
    set(gca, 'FontSize', 20); xlim([0 500])
    title(coordROI{gg}, 'Interpreter', 'none')
    
end
