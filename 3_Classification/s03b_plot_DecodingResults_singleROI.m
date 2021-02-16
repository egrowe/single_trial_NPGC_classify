%% PLOT ALL RESULTS TOGETHER

clear all; close all;

filePathCoords = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/3_Decoding/'];
filepath = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/3_Decoding/'];

for hh = 1:4
    
    if hh == 1 || hh == 4
        loadNameERP = ['AllEBResp_DiffVox_DecodeResults_SourcERP_Multifeature_' num2str(['L_Occ','L_FFA','L_InfPFC','L_SPPFC']) '.mat'];
        possCoords = {'L_Occ','L_FFA','R_Occ','R_FFA'}
    elseif hh == 2 || hh == 3
        loadNameERP = ['AllEBResp_DiffVox_DecodeResults_SourcERP_Multifeature_' num2str(['R_Occ','R_FFA','R_InfPFC','R_SPPFC']) '.mat'];
        possCoords = {'L_Occ','L_FFA','R_Occ','R_FFA'}
    end
    
    if hh == 1 || hh == 2
    load(loadNameERP, 'decodability_mean', 'decodability_std', 'allDecode')
    allERPdata(hh,:) = decodability_mean; allERPdata_std(hh,:) = decodability_std;
    allERPdata_sem(hh,:) =  allERPdata_std(hh,:)/sqrt(length(allDecode));
    end
    
    for gg = 1:length(possCoords)
        
        thisCoord = possCoords{gg}
        
        if hh == 1
            loadnameNonPara = ['resultsDecode_0to600msCohPowerGC_L_Inf_PFC_and_' num2str(thisCoord) '_TopVoxel_FINAL.mat'];            
        elseif hh == 2
            loadnameNonPara = ['resultsDecode_0to600msCohPowerGC_R_Inf_PFC_and_' num2str(thisCoord) '_TopVoxel_FINAL.mat'];            
        elseif hh == 3
            loadnameNonPara = ['resultsDecode_0to600msCohPowerGC_L_SP_PFC_and_' num2str(thisCoord) '_TopVoxel_FINAL.mat'];            
        elseif hh == 4
            loadnameNonPara = ['resultsDecode_0to600msCohPowerGC_R_SP_PFC_and_' num2str(thisCoord) '_TopVoxel_FINAL.mat'];            
        end
        
        load([filepath loadnameNonPara])
        results_NP = results; clear results; 
        data_NP = data; clear data;
        
        %Decodability results
        allPFC_logpwr(hh,gg,:) = results_NP.multiLogPwr_NodeB;
        allPFC_logpwr_std(hh,gg) = results_NP.multiLogPwr_NodeB_std;
        allPFC_logpwr_sem(hh,gg) = allPFC_logpwr_std(hh,gg)/sqrt(length(results_NP.multiLogPwr_NodeB));
        alllogPwr(hh,gg,:) = results_NP.multiLogPwr_NodeA;
        alllogPwr_std(hh,gg) = results_NP.multiLogPwr_NodeA_std;
        alllogPwr_sem(hh,gg) = alllogPwr_std(hh,gg)/sqrt(length(results_NP.multiLogPwr_NodeA));
        allCoh(hh,gg) = results_NP.multiCoh_mean;
        allCoh_sem(hh,gg) = results_NP.multiCoh_std;
        allCoh_sem(hh,gg) = allCoh_sem(hh,gg)/sqrt(length(results_NP.multiCoh));
        allGC_toPFC(hh,gg) = results_NP.multiGC_NodeA_to_NodeB_mean;
        allGC_toPFC_std(hh,gg) = results_NP.multiGC_NodeA_to_NodeB_std;
        allGC_toPFC_sem(hh,gg) = allGC_toPFC_std(hh,gg)/sqrt(length(results_NP.multiGC_NodeA_to_NodeB));
        allGC_fromPFC(hh,gg) = results_NP.multiGC_NodeB_to_NodeA_mean;
        allGC_fromPFC_std(hh,gg) = results_NP.multiGC_NodeB_to_NodeA_std;
        allGC_fromPFC_sem(hh,gg) = allGC_fromPFC_std(hh,gg)/sqrt(length(results_NP.multiGC_NodeB_to_NodeA));
        
        %Raw data 
        allGC_vals_toPFC_faces(hh,gg) = mean(mean(data_NP.GC_NodeA_to_NodeB_faces));
        allGC_vals_toPFC_std_faces(hh,gg) = std(mean(data_NP.GC_NodeA_to_NodeB_faces));
        allGC_vals_toPFC_sem_faces(hh,gg) = allGC_vals_toPFC_std_faces(hh,gg) /(sqrt(length(mean(data_NP.GC_NodeA_to_NodeB_faces))));
        allGC_vals_fromPFC_faces(hh,gg) = mean(mean(data_NP.GC_NodeB_to_NodeA_faces));
        allGC_vals_fromPFC_std_faces(hh,gg) = std(mean(data_NP.GC_NodeB_to_NodeA_faces));
        allGC_vals_fromPFC_sem_faces(hh,gg) = allGC_vals_fromPFC_std_faces(hh,gg)/(sqrt(length(mean(data_NP.GC_NodeA_to_NodeB_faces))));
        
        allGC_vals_toPFC_random(hh,gg) = mean(mean(data_NP.GC_NodeA_to_NodeB_random));
        allGC_vals_toPFC_std_random(hh,gg) = std(mean(data_NP.GC_NodeA_to_NodeB_random));
        allGC_vals_toPFC_sem_random(hh,gg) = allGC_vals_toPFC_std_random(hh,gg)/(sqrt(length(mean(data_NP.GC_NodeA_to_NodeB_faces))));
        allGC_vals_fromPFC_random(hh,gg) = mean(mean(data_NP.GC_NodeB_to_NodeA_random));
        allGC_vals_fromPFC_std_random(hh,gg) = std(mean(data_NP.GC_NodeB_to_NodeA_random));
        allGC_vals_fromPFC_sem_random(hh,gg) = allGC_vals_fromPFC_std_random(hh,gg)/(sqrt(length(mean(data_NP.GC_NodeA_to_NodeB_faces))));
    
         %allPwr_data_faces(hh,gg) = 
    end
    
    meanPFC_pwr(hh,:) = (allPFC_logpwr(hh,:));
    meanPFC_pwr_std(hh) = mean(allPFC_logpwr_std(hh,:));
    meanPFC_pwr_sem(hh) = mean(allPFC_logpwr_sem(hh,:));
    
end

% 
% %Gather all data to plot
% allPlotData = [allERPdata, meanPFC_pwr', allPwr, allCoh, allGC_toPFC, allGC_fromPFC];
% % 
% % %Gather std data to plot
% allPlotData_std = [allERPdata_std, meanPFC_pwr_std', allPwr_std, allCoh_std, allGC_toPFC_std, allGC_fromPFC_std];

%GATHER ALL POWER DATA
alllogPwrs = [meanPFC_pwr(1:4,1:10); squeeze(alllogPwr(1,:,:))];
alllogPwrs = [alllogPwrs(5,:); alllogPwrs(6,:);alllogPwrs(3,:);alllogPwrs(1,:);alllogPwrs(7,:);...
   alllogPwrs(8,:);alllogPwrs(4,:);alllogPwrs(2,:)]

alllogPwrs_sem = [meanPFC_pwr_sem(1:4), squeeze(alllogPwr_sem(1,:,:))]';
alllogPwrs_sem = [alllogPwrs_sem(5,:); alllogPwrs_sem(6,:);alllogPwrs_sem(3,:);alllogPwrs_sem(1,:);alllogPwrs_sem(7,:);...
   alllogPwrs_sem(8,:);alllogPwrs_sem(4,:);alllogPwrs_sem(2,:)]';


% REARRANGE ERP data to fit Occ, FFA, Par, ITG, PFC order
allERP_inOrder = allERPdata(:)'; %[allERPdata(1:2,5), allERPdata(1:2,2),allERPdata(1:2,4),allERPdata(1:2,3),allERPdata(1:2,1)]
allERP_sem_inOrder = allERPdata_sem(:)'; %[allERPdata_sem(1:2,5), allERPdata_sem(1:2,2),allERPdata_sem(1:2,4),allERPdata_sem(1:2,3),allERPdata_sem(1:2,1)]

%REARRANGE POWER data to fit Occ, FFA, Par, ITG, PFC order
allLogPwr_inOrder = mean(alllogPwrs,2); %[alllogPwrs(1:2,5), alllogPwrs(1:2,2),alllogPwrs(1:2,4),alllogPwrs(1:2,3),alllogPwrs(1:2,1)]
allLogPwr_sem_inOrder = alllogPwr_sem; %[alllogPwrs_sem(1:2,5), alllogPwrs_sem(1:2,2),alllogPwrs_sem(1:2,4),alllogPwrs_sem(1:2,3),alllogPwrs_sem(1:2,1)]

%REARRANGE COHERENCE data to fit Occ, FFA, Par, ITG, PFC order
allCoh_inOrder = allCoh; %[allCoh(:,2),allCoh(:,4),allCoh(:,3),allCoh(:,1)]
allCoh_sem_inOrder = allCoh_sem; %[allCoh_sem(:,2),allCoh_sem(:,4),allCoh_sem(:,3),allCoh_sem(:,1)]

%%REARRANGE GC to PFC data to fit Occ, FFA, Par, ITG, PFC order
allGC_toPFC_inOrder = allGC_toPFC; %[allGC_toPFC(:,4),allGC_toPFC(:,1),allGC_toPFC(:,3),allGC_toPFC(:,2)]
allGC_toPFC_sem_inOrder = allGC_toPFC_sem;%  [allGC_toPFC_sem(:,4),allGC_toPFC_sem(:,1),allGC_toPFC_sem(:,3),allGC_toPFC_sem(:,2)]


%%REARRANGE GC from PFC data to fit Occ, FFA, Par, ITG, PFC order
allGC_fromPFC_inOrder = allGC_fromPFC; %[allGC_fromPFC(:,4),allGC_fromPFC(:,1),allGC_fromPFC(:,3),allGC_fromPFC(:,2)]
allGC_fromPFC_sem_inOrder = allGC_fromPFC_sem; % [allGC_fromPFC_sem(:,4),allGC_fromPFC_sem(:,1),allGC_fromPFC_sem(:,3),allGC_fromPFC_sem(:,2)]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT ERP and log POWER (in one figure) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % GATHER ALL DATA TO PLOT FOR ERP and PWR decodability
plot_ERP_dataALL = allERP_inOrder;
plot_sem_ERP_dataALL = allERP_sem_inOrder;
% 
plot_pwr_ALL = allLogPwr_inOrder;
plot_sem_pwr_ALL = alllogPwrs_sem;

%Setup the labels
XLabels_pwr = {'Occ','FG','Inf-PFC','SP-PFC','Occ','FG','Inf-PFC','SP-PFC'}
    
    
% PLOT THE ERP SINGLE-ROI DECODING
figure(11)
errorbar(plot_ERP_dataALL(1,:)',plot_sem_ERP_dataALL(1,:)','b-','LineWidth', 5); hold on;
ylabel('Classification Accuracy (%)','FontSize',20);
set(gca, 'XTick', 1:length(plot_ERP_dataALL), 'XTickLabel', XLabels_pwr, 'FontSize',20);
xlim([0 length(plot_ERP_dataALL)+1]); ylim([48 90]);
%title(['Decodability Waveforms - Face v Random Trials'], 'FontSize', 30);
plot(0:length(plot_ERP_dataALL)+1, ones(1,length(plot_ERP_dataALL)+2)*50, '-r', 'LineWidth', 2); hold on;
legend({'ERP'})


figure(22)

%PLOT THE POWER SINLE ROI-DECODING
subplot(2,1,1)
violinplot(alllogPwrs', XLabels_pwr)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT COHERENCE and POWER DECODING RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allCoh_plot = allCoh_inOrder;
allCoh_sem = allCoh_sem_inOrder;

%Setup the labels
XLabels = {'Occ','FG','Occ','FG'}


figure(22)

%PLOT THE POWER SINLE ROI-DECODING
subplot(2,1,1)
errorbar(plot_pwr_ALL,plot_sem_pwr_ALL,'k-', 'LineWidth', 5); hold on;
ylabel('Classification Accuracy (%)','FontSize',20);
set(gca, 'XTick', 1:length(plot_pwr_ALL), 'XTickLabel', XLabels_pwr, 'FontSize',20);
xlim([0 length(plot_pwr_ALL)+1]); ylim([48 90]);
%title(['Decodability Waveforms - Face v Random Trials'], 'FontSize', 30);
plot(0:length(plot_pwr_ALL)+1, ones(1,length(plot_pwr_ALL)+2)*50, '-r', 'LineWidth', 2); hold on;
legend({'Pwr'})

%PLOT THE COHERNECE SINGLE ROI-DECODING
subplot(2,1,2)
errorbar(allCoh_plot(1,:),allCoh_sem(1,:),'b-','LineWidth', 5); hold on;
errorbar(allCoh_plot(2,:),allCoh_sem(2,:), 'k-','LineWidth', 5); hold on;
errorbar(allCoh_plot(3,:),allCoh_sem(3,:),'k:','LineWidth', 5); hold on;
errorbar(allCoh_plot(4,:),allCoh_sem(4,:), 'b:','LineWidth', 5); hold on;
ylabel('Classification Accuracy (%)','FontSize',20);
set(gca, 'XTick', 1:length(allCoh_plot), 'XTickLabel', XLabels, 'FontSize',20);
xlim([0 length(allCoh_plot)+1]); ylim([45 65]);
%title(['Decodability Coherence - Face v Random Trials'], 'FontSize', 20);
plot(0:length(allCoh_plot)+1, ones(1,length(allCoh_plot)+2)*50, '-r', 'LineWidth', 2); hold on;
legend({'Left Inf PFC and ROI','Right Inf PFC and ROI','Right SP PFC and ROI','Left SP PFC and ROI'}, 'FontSize',15)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT Granger Causality DECODING RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_GC_to_PFC_ALL = allGC_toPFC_inOrder
plot_sem_GC_to_PFC_ALL = allGC_toPFC_sem_inOrder

plot_GC_from_PFC_ALL = allGC_fromPFC_inOrder;
plot_sem_GC_from_PFC_ALL = allGC_fromPFC_sem_inOrder;

XLabels = {'L Occ','L FG','R Occ','R FG'}

figure(212)
subplot(2,1,2)
errorbar(plot_GC_from_PFC_ALL(1,:),plot_sem_GC_from_PFC_ALL(1,:),'LineWidth', 5); hold on;
errorbar(plot_GC_from_PFC_ALL(2,:),plot_sem_GC_from_PFC_ALL(2,:), 'LineWidth', 5); hold on;
errorbar(plot_GC_from_PFC_ALL(3,:),plot_sem_GC_from_PFC_ALL(3,:),'LineWidth', 5); hold on;
errorbar(plot_GC_from_PFC_ALL(4,:),plot_sem_GC_from_PFC_ALL(4,:), 'LineWidth', 5); hold on;
ylabel('Classification Accuracy (%)','FontSize',20);
set(gca, 'XTick', 1:length(plot_GC_from_PFC_ALL), 'XTickLabel', XLabels, 'FontSize',20);
xlim([0 length(plot_GC_from_PFC_ALL)+1]); ylim([42 70]);
%title(['Decodability Granger Causality - Face v Random Trials'], 'FontSize', 20);
plot(0:length(plot_GC_from_PFC_ALL)+1, ones(1,length(plot_GC_from_PFC_ALL)+2)*50, '-r', 'LineWidth', 2); hold on;
legend({'Left Inf PFC to ROI','Right Inf PFC to ROI','Right SP PFC to ROI','Left SP PFC to ROI'}, 'FontSize',15)

subplot(2,1,1)
errorbar(plot_GC_to_PFC_ALL(1,:),plot_sem_GC_to_PFC_ALL(1,:),'LineWidth', 5); hold on;
errorbar(plot_GC_to_PFC_ALL(2,:),plot_sem_GC_to_PFC_ALL(2,:), 'LineWidth', 5); hold on;
errorbar(plot_GC_to_PFC_ALL(3,:),plot_sem_GC_to_PFC_ALL(3,:),'LineWidth', 5); hold on;
errorbar(plot_GC_to_PFC_ALL(4,:),plot_sem_GC_to_PFC_ALL(4,:), 'LineWidth', 5); hold on;
ylabel('Classification Accuracy (%)','FontSize',20);
set(gca, 'XTick', 1:length(plot_GC_to_PFC_ALL), 'XTickLabel', XLabels, 'FontSize',20);
xlim([0 length(plot_GC_to_PFC_ALL)+1]); ylim([42 70]);
%title(['Decodability Granger Causality - Face v Random Trials'], 'FontSize', 20);
plot(0:length(plot_GC_to_PFC_ALL)+1, ones(1,length(plot_GC_to_PFC_ALL)+2)*50, '-r', 'LineWidth', 2); hold on;
legend({'ROI to Left Inf PFC','ROI to Right Inf PFC','ROI to Right SP PFC','ROI to Left SP PFC'}, 'FontSize',15)

