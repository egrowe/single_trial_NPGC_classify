%% PLOT ALL RESULTS TOGETHER

clear all; close all;
addpath(genpath('/Users/egrow1/Desktop/Shafto_GitHub_Scripts/3_Decoding/RainCloudPlots-master'));
addpath('/Users/egrow1/Desktop/Shafto_GitHub_Scripts/i_GC_Simulation/mseb');

%ROI combinations to plot decoding accuracy
possROICombs = {'ALL_PFC_ALL_FG_OCC','ALL_PFC_ALL_FG',...
   'ALL_PFC_ALL_OCC'}; %'main' PFC coords
maxFreq = 55;


%possROICombs = {'ALL_FG_ALL_OCC'}; %'main' PFC coords

for gg = 1:length(possROICombs)
    
    useThisComb = possROICombs{gg}
    

       loadname = ['results_noF_Decode_CohPowerGC_Coords_' num2str(maxFreq) 'Hz_MultiROI_' ...
           num2str(useThisComb) '_5Tapsmo_TopVoxel_FINAL.mat']; %savename
            
    
    load(loadname)
    results_NP = results; clear results;
    data_NP = data; clear data;
    
    %Decodability results
    allGC_toPFC(gg,:) = results_NP.multiGC_NodeA_to_NodeB;
    allGC_fromPFC(gg,:) = results_NP.multiGC_NodeB_to_NodeA;
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT Granger Causality DECODING RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT USING RAINCLOUD
clear newData
% get nice colours from colorbrewer
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');

cl(1,:) = cb(5, :);
cl(2,:) = cb(6, :);
cl(3,:) = cb(8, :);

% WHICH DATA TO PLOT?
for bb = 1:size(allGC_toPFC,1)
    newData{bb,:} = allGC_toPFC(bb,:)'
end

% %Plot
% make figure
figure(1)
subplot(1,2,1)
h   = rm_raincloud(newData, cl(1,:));
set(gca, 'FontSize',20)

% PLOT BACK FROM PFC TO SENSORY ROIS
% WHICH DATA TO PLOT?
for bb = 1:size(allGC_fromPFC,1)
    newData{bb,:} = allGC_fromPFC(bb,:)'
end

% %Plot
% make figure
figure(1)
subplot(1,2,2)
h   = rm_raincloud(newData, cl(1,:));
set(gca, 'FontSize',20)
title('GC decoding results')
