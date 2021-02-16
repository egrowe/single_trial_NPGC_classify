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
    allPFC_logpwr(gg,:) = results_NP.multiLogPwr_NodeB;
    
    allPFC_data_logPwr_faces = data_NP.all_PFClogpwr_faces;
    allOther_data_logPwr_faces = data_NP.all_ROIlogpwr_faces;
    allPFC_data_logPwr_random = data_NP.all_PFClogpwr_random;
    allOther_data_logPwr_random = data_NP.all_ROIlogpwr_random ;
    
    
    alllogPwr(gg,:) = results_NP.multiLogPwr_NodeA;
    allCoh(gg,:) = results_NP.multiCoh;
    allGC_toPFC(gg,:) = results_NP.multiGC_NodeA_to_NodeB;
    allGC_fromPFC(gg,:) = results_NP.multiGC_NodeB_to_NodeA;
    
    if gg == 1
        %% PLOT DATA FOR POWER
        allPwr_faces = [allPFC_data_logPwr_faces; allOther_data_logPwr_faces];
        allPwr_random = [allPFC_data_logPwr_random; allOther_data_logPwr_random];
        
        for vv = 1:size(allPwr_faces,1)
           
            minVal = min(allPwr_faces(:)); freqWindows = linspace(0,251,127);
            %Face data
            LogPwr_faces = squeeze(allPwr_faces(vv,:,:));
            mLogPwr_faces = mean(LogPwr_faces,2);
            std_Logpwr_faces = std(LogPwr_faces')
            sem_Logpwr_faces = std_Logpwr_faces/(sqrt(length(LogPwr_faces))); %SEM of log transformed power
            
            %Random data
            LogPwr_random = squeeze(allPwr_random(vv,:,:));
            mLogPwr_random = mean(LogPwr_random,2);
            std_Logpwr_random = std(LogPwr_random')
            sem_Logpwr_random = std_Logpwr_random/(sqrt(length(LogPwr_random))); %SEM of log transformed power
            
            %Plot
            figure(88)
            subplot(4,2,vv)
            errorbar(mLogPwr_faces, sem_Logpwr_faces, 'LineWidth', 2); hold on;
            errorbar(mLogPwr_random,sem_Logpwr_random, 'LineWidth', 2); hold on;
            ylim([minVal 0.2]); xlim([1 size(LogPwr_faces,1)]); % ylim([-3 0]); %ylabel('Coherence','FontSize',12);
           if vv == 7
                       set(gca, 'XTick', [1,round(length(freqWindows)/2),length(freqWindows)], 'XTickLabel', ...
                {ceil(freqWindows(1)),ceil(freqWindows(ceil(length(freqWindows)/2))),ceil(freqWindows(end))}); 
            set(gca, 'FontSize',20); ylabel('log Power'); xlabel('Frequency (Hz)');
           else
            set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]);
           end
          
           allPwr = [mLogPwr_faces'; mLogPwr_random'];
           allPwr_sem = [sem_Logpwr_faces; sem_Logpwr_random];
           figure(818181)
                 subplot(4,2,vv)
        mseb((repmat(1:length(allPwr),2,1)),allPwr,allPwr_sem,[],1); box on;
                    ylim([-3.5 -0.5]); xlim([1 size(LogPwr_faces,1)]); % ylim([-3 0]); %ylabel('Coherence','FontSize',12);      
        if vv == 7
        set(gca, 'XTick', [1,round(length(freqWindows)/2),length(freqWindows)], 'XTickLabel', ...
                {ceil(freqWindows(1)),ceil(freqWindows(ceil(length(freqWindows)/2))),ceil(freqWindows(end))}); 
            set(gca, 'FontSize',20); ylabel('log Power'); xlabel('Frequency (Hz)');
           else
            set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]);
           end
        end
    end
    
end

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT Coherence DECODING RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT COHERNECE BETWEEN PFC and SENSORY ROIs
clear newData
% WHICH DATA TO PLOT?
for bb = 1:size(allCoh,1)
    newData{bb,:} = allCoh(bb,:)'
end

% %Plot
% make figure
figure(2)
h   = rm_raincloud(newData, cl(1,:));
set(gca, 'FontSize',20)
title('Coherence decoding')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT COMBINED Log Power DECODING RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gather data for "ALL"  pwr
allComb_logPwr = results_NP.multiLogPwr_NodeA_and_B;

% WHICH DATA TO PLOT?
clear newData
newData{1,:} = allComb_logPwr'
newData{2,:} = allPFC_logpwr(1,:)'
newData{3,:} = alllogPwr(1,:)'


% %Plot
% make figure
figure(3)
h   = rm_raincloud(newData, cl(1,:));
set(gca, 'FontSize',20)
title('Log power decoding results')





%% EXTRACT INDIVIDUAL ROI POWER DATA
clear allPFC_logpwr alllogPwr
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');

cl(1, :) = cb(4, :);
cl(2, :) = cb(1, :);
% SELECT THE NODES TO EXTRACT
useROIs = {'L_PFC_and_L_FFA','R_PFC_and_R_FFA','L_Inf_PFC_and_L_Occ','R_Inf_PFC_and_R_Occ'}

for bb = 1:length(useROIs)
    
    usethisROI = useROIs{bb}
    
    loadThis = ['results_NoF_Decode_0to500ms_CohPowerGC_' num2str(usethisROI) ...
        '_' num2str(maxFreq) 'Hz_5Tapsmo_TopVoxel_FINAL.mat']; %savename
    
    load(loadThis, 'results','data')
    results_NP = results; clear results;
    data_NP = data; clear data;
    
    %Decodability results
    allPFC_lpwr(bb,:) = results_NP.multiLogPwr_NodeB
    all_lPwr(bb,:) = results_NP.multiLogPwr_NodeA
    
    clear results_NP data_NP
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT log power DECODING RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Distribute into columns for LEFT and RIGHT hemispheres
for gg = 1:size(allPFC_lpwr,1)
    
    if gg == 1
        newData{1,1} = allPFC_lpwr(gg,:)';
        newData{3,1} = all_lPwr(gg,:)';
    elseif gg == 2
        newData{1,2} = allPFC_lpwr(gg,:)';
        newData{3,2} = all_lPwr(gg,:)';
    elseif gg == 3
        newData{2,1} = allPFC_lpwr(gg,:)';
        newData{4,1} = all_lPwr(gg,:)';
    elseif gg == 4
        newData{2,2} = allPFC_lpwr(gg,:)';
        newData{4,2} = all_lPwr(gg,:)';
    end
end

%% PLOT MULTI PLOTS (connected)
% read into cell array of the appropriate dimensions
% make figure
figure(4)
h   = rm_raincloud(newData, cl(1:2,:));
set(gca, 'FontSize',20)
title('Log power (individual) decoding results')
