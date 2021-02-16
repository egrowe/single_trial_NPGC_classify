% %% Extract individual ERPs and plot over time
clear all; close all;

% Add the libSVM decoding library to the path
addpath(genpath('/Users/egrow1/Desktop/NEW_Redo_OriginalChans_SHAFTO/decode_libSVM'))
addpath(genpath('/Users/egrow1/Desktop/Shafto_GitHub_Scripts/3_Decoding/RainCloudPlots-master'))

filepath = ['/Users/egrow1/Desktop/ALL_SW_noF_0_to_500ms_EBR/']; %set file path
filepathData = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/1_EEG_PreProcessing/']; %coordinates

% Where are we extractig the data from?
hemisphere = 'Both'; % which hemisphere? Either 'Left' or 'Right', 'Both' (i.e all) or 'CombHemi' (combined each ROI bilaterally)

%Assign other parameters
zscore = 1; %1 = on, 0 = off
timeWindow = [0 500]; %examine this time window
timeIdx = linspace(0,500,251); % find index for time points in data
useWind = [find(timeIdx==timeWindow(1)):find(timeIdx==timeWindow(2))]; %assign idxs by time window
visualise = 0; %to show data visualisation =1, =0 for not
count = 1; %set up counter to plot (if visualise is turned on)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine the coordinates to extract ROI data from
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assign variables depending on which hemisphere being examined
if strcmp(hemisphere, 'Left')
    %Name the left/right coordinates to examine and labels and savename
    coordLabels = {'L_PFC','L_Inf_PFC','L_FFA','L_Occ'};
     savename = ['Decode_SourcERP_Multifeature_Hemispheres=' num2str(hemisphere) '.mat'];
elseif strcmp(hemisphere, 'Right')
    coordLabels = {'R_PFC','R_Inf_PFC','R_FFA','R_Occ'};
     savename = ['Decode_SourcERP_Multifeature_Hemispheres=' num2str(hemisphere) '.mat'];
elseif strcmp(hemisphere, 'Both')
    coordLabels = {'L_PFC','L_Inf_PFC','L_FFA','L_Occ',...
        'R_PFC','R_Inf_PFC','R_FFA','R_Occ'};
     savename = ['Decode_SourcERP_Multifeature_Hemispheres=' num2str(hemisphere) '.mat'];
elseif strcmp(hemisphere, 'CombHemi')
    coordLabels = {'L_PFC','R_PFC', 'L_Inf_PFC', 'R_Inf_PFC', ...
        'L_FFA', 'R_FFA', 'L_Occ', 'R_Occ'}; 
     savename = ['Decode_SourcERP_Multifeature_Hemispheres=' num2str(hemisphere) '.mat'];
end

% Organise coordinates if decoding bilaterally (i.e. 'CombHemi')
cc = 1; ff = 1;
for gg = 1:length(coordLabels)
    
    load([filepathData 'Coord_NoF_FINAL_' num2str(coordLabels{gg}) '_by_TstatOverTime_0_to_500ms'])

    if  strcmp(hemisphere, 'CombHemi')
        if cc <= ff
        useCoords{cc,1} = sortedCoords(1,:); cc=cc+1;
        else
        useCoords{ff,2} = sortedCoords(1,:); ff=ff+1;
        end
    else  
    useCoords{gg,:} = sortedCoords(1,:)
    end

end

%% EXTRACT DATA FROM EACH SOURCE COORDINATE
for gg = 1:length(useCoords)
    
    coordList = useCoords{gg}
    
    %% Do it for face trials
    for voxel = 1:size(coordList,1)
        thisVoxel = coordList;
        filename = ['NEW_Source_Waveform_0to500ms_at_' num2str(thisVoxel) '_for_ALL_trials_mSae_faces_aEBfspmeeg_7MJN_0to500ms_trial_205.mat'];
        load ([filepath filename])
        GM_faces(voxel,:,:) = waveform;%(:,useWind); % time x trial x feature
    end
    
    %% Do it for RANDOM trials
    for voxel = 1:size(coordList,1)
        thisVoxel = coordList;
        filename = ['NEW_Source_Waveform_0to500ms_at_' num2str(thisVoxel) '_for_ALL_trials_mSae_random_aEBfspmeeg_7MJN_0to500ms_trial_354.mat'];
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
                hemisphere filepathResp zscore useWind visualise coordLabels savename count saveFigDecode saveFigWaves
    
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

%% PLOT USING RAINCLOUD
f4 = figure;
% get nice colours from colorbrewer
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');

cl(1, :) = cb(4, :);
cl(2, :) = cb(1, :);

fig_position = [200 200 600 400]; 

% WHICH DATA TO PLOT?
%  
% for bb = 1:size(allDecode,1)
%     newData{bb,:} = allDecode(bb,:)'
% end

%Distribute into columns for LEFT and RIGHT hemispheres
cc = 1; ff = 1; ss = 1;
for gg = 1:length(coordLabels)
    if  strcmp(hemisphere, 'Both')
        if cc <= ff
        newData{cc,1} = allDecode(ss,:)'; cc=cc+1; ss=ss+1;
        else
        newData{ff,2} = allDecode(ss,:)'; ff=ff+1; ss=ss+1;
        end
    else  
    useCoords{gg,:} = sortedCoords(1,:)
    end
end


%% PLOT MULTI PLOTS (connected)
% read into cell array of the appropriate dimensions
% make figure
h   = rm_raincloud(newData, cl);
set(gca, 'FontSize',20)
title(['Source waveform decoding: All ROIs (individually)']);


%Save these results
save(savename)
