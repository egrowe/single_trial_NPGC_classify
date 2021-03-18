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
filepath = ['/Users/egrow1/Desktop/ALL_SourceWaves_8196_EBRemoved/']; %set file path
filepathB = ['/Users/egrow1/Desktop/Shafto_Analysis_Scripts/Single-Trial/EEG_and_Source_PreProcessing/']
%filepathRand = ['/Users/egrow1/Desktop/ALL_new_SingleBins_359RandomTrials/'];
%filepathFace = ['/Users/egrow1/Desktop/ALL_NEW_Sourcewaves_NewTrials_8196Voxels/'];
filepathCoords = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/1_EEG_PreProcessing/'];
%filepathResp = ['/Users/egrow1/Desktop/ALL_new_faceResp_SourceWaves/'];
listCoords = ['Coordinates_FINAL_UPDATED_8196_voxels_byROI.mat'];
load([filepathCoords listCoords])
coordList = L_Occ;
nameROI = 'L_Occ';

%Define other parameters
zscore = 1; %1 = on, 0 = offclea
timeWindow = [0 600]; %examine this time window
timeIdx = linspace(0,600,301); % find index for time points in data
useWind = [find(timeIdx==timeWindow(1)):find(timeIdx==timeWindow(2))]; %assign idxs by time window
visualise = 1; %to show data visualisation =1, =0 for not

%Set the savename
savename = ['CoordSortedEBRem_' num2str(nameROI) '_by_TstatOverTime_' ...
    num2str(timeWindow(1)) '_to_' num2str(timeWindow(2)) 'ms.mat'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract waveforms for face trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for voxel = 1:length(coordList)
    thisVoxel = coordList(voxel,:);
    filename = ['SourceWaveform_at_' num2str(thisVoxel) '_allTrials_mSae_faces_aEBfspmeeg_P0291_0to500ms.mat'];
    load ([filepath filename])
    GM_faces(voxel,:,:) = waveform;%(:,useWind); % time x trial x feature
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract waveforms for RANDOM trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for voxel = 1:length(coordList)
    thisVoxel = coordList(voxel,:);
    filename = ['SourceWaveform_at_' num2str(thisVoxel) '_allTrials_mSae_random_aEBfspmeeg_P0291_0to500ms.mat'];
    load ([filepath filename])
    GM_random(voxel,:,:) = waveform;%(:,useWind); % time x trial x feature
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gather all data within time window AND VISUALISE (if required)
% GATHER DATA FOR DECODING
GM_faces = GM_faces(:,:,useWind); %restrict to specified time window (if applicable)
GM_random = GM_random(:,:,useWind); %restrict to specified time window (if applicable)

%% Plot the timecourse
for voxel = 1:length(coordList)
    %plot timecourse of tstat
    for ii = 1:size(GM_faces,3)
        
        thistimePoint_F = squeeze(GM_faces(voxel,:,ii));
        thisTimePoint_NF = squeeze(GM_random(voxel,:,ii));
        
        [h,p,i,stat] = ttest2(thistimePoint_F,thisTimePoint_NF);
        
        tstat(voxel,ii) = stat.tstat;
    end
end

% %Plot before sorting
figure(4);
tstat = abs(tstat)
imagesc(tstat, [0 8.5]); hold on; colorbar;
xlabel('Time (ms)'); ylabel('Voxel #');
set(gca,'XTick',[1, 151, 301], 'XTickLabel', {'0','300','600'});
set(gca,'FontSize', 20)
%
% %Sort by maximum t-stat value over time
maxEachRow = max(tstat')
[valMax,y] = sort(maxEachRow, 'descend')
sortedTstat = tstat(y,:);
%
% % %Plot after sorting
figure(5)
imagesc(sortedTstat, [0 8.5]); hold on; colorbar;
xlabel('Time (ms)'); ylabel('Voxel #');
set(gca,'XTick',[1, 151, 301], 'XTickLabel', {'0','300','600'});
set(gca,'FontSize', 20)
title(['T state over time, sorted, for all voxels in ROI:' num2str(nameROI)], 'FontSize', 15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SORT THE VOXEL COORDS BY HIGHEST T-STAT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sort coordinates by tstat over time
sortedCoords = coordList(y,:)
topCoord = sortedCoords(1,:)

save(savename, 'sortedTstat','tstat','sortedCoords')