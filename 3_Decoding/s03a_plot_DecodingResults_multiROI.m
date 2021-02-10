%% PLOT ALL RESULTS TOGETHER

clear all; close all;

hh = 1;

possMainCoords = {'_ALL_PFC_ALL_FG_OCC_','_ALL_PFC_ALL_FG_',...
    '_ALL_PFC_ALL_OCC_'}; %'main' PFC coords

for gg = 1:length(possMainCoords)
    
    useMainCoord = possMainCoords{gg}

        loadnameNonPara = ['results_noF_Decode_CohPowerGC_CoordsMultiROI' ...
            num2str(useMainCoord) 'TopVoxel_FINAL.mat']; %savename
        
        load(loadnameNonPara)
        results_NP = results; clear results;
        data_NP = data; clear data;
        
        %Decodability results
        allPFC_logpwr(gg,:) = results_NP.multiLogPwr_NodeB;
        alllogPwr(gg,:) = results_NP.multiLogPwr_NodeA;
        allCoh(gg) = results_NP.multiCoh_mean;
        allGC_toPFC(gg) = results_NP.multiGC_NodeA_to_NodeB_mean;
        allGC_fromPFC(gg) = results_NP.multiGC_NodeB_to_NodeA_mean;
      
end

meanPFC_pwr = squeeze(mean(allPFC_logpwr,2));

%% GATHER ALL RESULTS TO PLOT
%
% %Gather all data to plot
allPlotData = [meanPFC_pwr, alllogPwr, allCoh', allGC_toPFC', allGC_fromPFC'];
%
% %Gather std data to plot

%GATHER ALL POWER DATA
alllogPwrs = [meanPFC_pwr, alllogPwr];

alllogPwrs = [alllogPwrs(2,:); alllogPwrs(1,:); alllogPwrs(3,:); alllogPwrs(4,:)]

%REARRANGE POWER data to fit Occ, FFA, Par, ITG, PFC order
allLogPwr_inOrder = alllogPwrs; % [alllogPwrs(:,5), alllogPwrs(:,2),alllogPwrs(:,4),alllogPwrs(:,3),alllogPwrs(:,1)]

%REARRANGE COHERENCE data to fit Occ, FFA, Par, ITG, PFC order
allCoh_inOrder =[allCoh(2,:),allCoh(:,4),allCoh(:,3),allCoh(:,1)]

%%REARRANGE GC to PFC data to fit Occ, FFA, Par, ITG, PFC order
allGC_toPFC_inOrder = allGC_toPFC; % [allGC_toPFC(:,4),allGC_toPFC(:,1),allGC_toPFC(:,3),allGC_toPFC(:,2)]

%%REARRANGE GC from PFC data to fit Occ, FFA, Par, ITG, PFC order
allGC_fromPFC_inOrder = allGC_fromPFC; %[allGC_fromPFC(:,4),allGC_fromPFC(:,1),allGC_fromPFC(:,3),allGC_fromPFC(:,2)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT ERP and log POWER (in one figure) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GATHER ALL DATA TO PLOT FOR ERP and PWR decodability
plot_ERP_dataALL = allERP_inOrder;
plot_sem_ERP_dataALL = allERP_sem_inOrder;

plot_pwr_ALL = allLogPwr_inOrder;
plot_sem_pwr_ALL = allLogPwr_sem_inOrder;

%Setup the labels
XLabels = {'All ROIs','All PFC','All Occ','All FG'}
figure(21)
subplot(2,1,1)
violinplot(alllogPwrs', XLabels)


%errorbar(plot_ERP_dataALL(1,:)',plot_sem_ERP_dataALL(1,:)','b-','LineWidth', 8); hold on;
errorbar(plot_pwr_ALL(1,:)',plot_sem_pwr_ALL(1,:)','k-','LineWidth', 8); hold on;
ylabel('Classification Accuracy (%)','FontSize',20);
set(gca, 'XTick', 1:length(plot_pwr_ALL), 'XTickLabel', XLabels, 'FontSize',20);
xlim([0 length(plot_pwr_ALL)+1]); ylim([48 90]);
%title(['Decodability Waveforms and log Power - Face v Random Trials'], 'FontSize', 20);
plot(0:length(plot_pwr_ALL)+1, ones(1,length(plot_pwr_ALL)+2)*50, '-r', 'LineWidth', 2); hold on;
legend({'log Power'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT COHERENCE DECODING RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allCoh_plot = allCoh;
allCoh_sem = allCoh_sem;

%Setup the labels
XLabels = {'ALL ROIs','All Occ','All FG'}

subplot(2,1,2)
errorbar(allCoh_plot(1,:),allCoh_sem(1,:),'b-','LineWidth', 6); hold on;
ylabel('Classification Accuracy (%)','FontSize',20);
set(gca, 'XTick', 1:length(allCoh_plot), 'XTickLabel', XLabels, 'FontSize',20);
xlim([0 length(allCoh_plot)+1]); ylim([45 60]);
%title(['Decodability Coherence - Face v Random Trials'], 'FontSize', 20);
plot(0:length(allCoh_plot)+1, ones(1,length(allCoh_plot)+2)*50, '-r', 'LineWidth', 2); hold on;
legend({'with multi PFCs'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT Granger Causality DECODING RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_GC_to_PFC_ALL = allGC_toPFC_inOrder
plot_sem_GC_to_PFC_ALL = allGC_toPFC_sem_inOrder

plot_GC_from_PFC_ALL = allGC_fromPFC_inOrder;
plot_sem_GC_from_PFC_ALL = allGC_fromPFC_sem_inOrder;

XLabels = {'ALL ROIs','All Occ','All FG'}

figure(212)
subplot(2,1,2)
errorbar(plot_GC_from_PFC_ALL(1,:),plot_sem_GC_from_PFC_ALL(1,:),'LineWidth', 5); hold on;
ylabel('Classification Accuracy (%)','FontSize',20);
set(gca, 'XTick', 1:length(plot_GC_from_PFC_ALL), 'XTickLabel', XLabels, 'FontSize',20);
xlim([0 length(plot_GC_from_PFC_ALL)+1]); ylim([45 70]);
%title(['Decodability Granger Causality - Face v Random Trials'], 'FontSize', 30);
plot(0:length(plot_GC_from_PFC_ALL)+1, ones(1,length(plot_GC_from_PFC_ALL)+2)*50, '-r', 'LineWidth', 2); hold on;
legend({'PFC to ROI'})

subplot(2,1,1)
errorbar(plot_GC_to_PFC_ALL(1,:),plot_sem_GC_to_PFC_ALL(1,:),'LineWidth', 5); hold on;
ylabel('Classification Accuracy (%)','FontSize',20);
set(gca, 'XTick', 1:length(plot_GC_to_PFC_ALL), 'XTickLabel', XLabels, 'FontSize',20);
xlim([0 length(plot_GC_to_PFC_ALL)+1]); ylim([45 70]);
%title(['Decodability Granger Causality - Face v Random Trials'], 'FontSize', 30);
plot(0:length(plot_GC_to_PFC_ALL)+1, ones(1,length(plot_GC_to_PFC_ALL)+2)*50, '-r', 'LineWidth', 2); hold on;
legend({'ROI to PFC'})
