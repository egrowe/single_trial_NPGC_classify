% %% Extract individual ERPs and plot over time
clear all; close all;

% Add the libSVM decoding library to the path
addpath(genpath('/Users/egrow1/Desktop/NEW_Redo_OriginalChans_SHAFTO/decode_libSVM'))

filepath = ['/Users/egrow1/Desktop/ALL_SW_noF_0_to_500ms_EBR/']; %set file path
filepathData = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/1_EEG_PreProcessing/']; %coordinates
savename = ['Decode_SourcERP_Multifeature.mat'];

%Assign other parameters
zscore = 1; %1 = on, 0 = off
timeWindow = [0 500]; %examine this time window
timeIdx = linspace(0,500,251); % find index for time points in data
useWind = [find(timeIdx==timeWindow(1)):find(timeIdx==timeWindow(2))]; %assign idxs by time window
visualise = 0; %to show data visualisation =1, =0 for not
count = 1; %set up counter to plot (if visualise is turned on)
nCombs = 3; %3 combinations of multi ROI decoding

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine the coordinates to extract ROI data from
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nn = 1:nCombs
 
    % WHICH ROIs are we decoding?
    if nn == 1
        coordLabels = {'L_PFC','L_Inf_PFC','L_FFA','L_Occ',...
            'R_PFC','R_Inf_PFC','R_FFA','R_Occ'};
    elseif nn == 2
        coordLabels = {'L_FFA','R_FFA'};
    elseif nn == 3
        coordLabels = {'L_Occ','R_Occ'};
    end
    
    %Extract coords for the 'top voxel' from these ROIs
    for gg = 1:length(coordLabels)
        load([filepathData 'Coord_NoF_FINAL_' num2str(coordLabels{gg}) '_by_TstatOverTime_0_to_500ms'])
        useCoords{gg,:} = sortedCoords(1,:)
    end
    
    %% EXTRACT DATA FROM EACH SOURCE COORDINATE
    for gg = 1:length(useCoords)
        
        coordList = useCoords{gg}
        
        %% Do it for face trials
        for voxel = 1:size(coordList,1)
            thisVoxel = coordList;
            filename = ['NEW_Source_Waveform_0to500ms_at_' num2str(thisVoxel) '_for_ALL_trials_mSae_faces_aEBfspmeeg_7MJN_0to500ms_trial_205.mat'];
            load ([filepath filename])
            GM_faces(gg,:,:) = waveform;%(:,useWind); % time x trial x feature
        end
        
        %% Do it for RANDOM trials
        for voxel = 1:size(coordList,1)
            thisVoxel = coordList;
            filename = ['NEW_Source_Waveform_0to500ms_at_' num2str(thisVoxel) '_for_ALL_trials_mSae_random_aEBfspmeeg_7MJN_0to500ms_trial_354.mat'];
            load ([filepath filename])
            GM_random(gg,:,:) = waveform;%(:,useWind); % time x trial x feature
        end
        
        
    end
    
    %% Gather all data within time window AND VISUALISE (if turned on)
    % GATHER DATA FOR DECODING
    data_faces = reshape(GM_faces, [size(GM_faces,1)*size(GM_faces,2), size(GM_faces,3)])'; %form features x trials
    data_random = reshape(GM_random, [size(GM_random,1)*size(GM_random,2), size(GM_random,3)])'; %form features x trials
    
    
    %% Decode using libSVM
    nReps = 10;
    costRange = -5;
    
    [decodability, cv_acc, w] = decode_multiFeature_libSVM(data_faces, data_random, nReps, costRange, zscore)
    
    %% Plot these results
    %Find mean and std of decodability
    decodability_all(nn,:) = decodability;

    clear decodability cv_acc w data_faces data_random
    
end

%% PLOT USING RAINCLOUD
f4 = figure;
% get nice colours from colorbrewer
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');

cl(1,:) = cb(5, :);
cl(2,:) = cb(6, :);
cl(3,:) = cb(8, :);

fig_position = [200 200 600 400];

% WHICH DATA TO PLOT?
for bb = 1:size(decodability_all,1)
    newData{bb,:} = decodability_all(bb,:)'
end

% %Plot
% make figure
% % define new colour
% new_cl_2 = cl(2,:);
% new_cl_3 = cl(3,:);
% 
% % change one subset to new colour and alter dot size
% h.p{2, 2}.FaceColor         = new_cl_2;
% h.s{2, 2}.MarkerFaceColor   = new_cl_2;
% h.m(2, 1).MarkerEdgeColor   = 'none';
% h.m(2, 1).MarkerFaceColor   = new_cl_2;
% h.s{2, 2}.SizeData          = 300;



h   = rm_raincloud(newData, cl(1,:));
set(gca, 'FontSize',20)
title(['Source waveform decoding: All ROIs (combined)']);



%
%Save these results
save(savename)