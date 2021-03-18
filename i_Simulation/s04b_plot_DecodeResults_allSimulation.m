%% PLOT ALL SIMULATION RESULTS

clear all; close all;

filePathCoords = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/3_Decoding/'];

%A matrix connectivity
poss_a = -0.5: 0.05: 0.5; %weight for 'a' (= x self conn X)
poss_c = -0.5: 0.05: 0.5; %weight for 'c' (= x self conn X in Y data)
poss_d = -0.5: 0.05: 0.5;  %weight for 'd' (= y self conn Y)

count = 0; nReps = 10;

for hh = 1:length(poss_a)
    
    for  ff = 1:length(poss_c)
        
        for vv = 1:length(poss_d)
            count = count + 1;
            
            %Model parameters
            this_a = poss_a(ss); % use this weight for 'a' (= x self conn X)
            this_c = poss_c(ff);  % use this weight for 'c' (= x self conn X in Y data)
            this_d = poss_d(vv); % use this weight for 'd' (= y self conn Y)
            
            loadname = ['All_Decode_Synthetic_SourceERP_multi_a=' num2str(this_a) '_c=' num2str(this_c) '_d=' num2str(this_d) '.mat']
            load([filePathCoords loadname], 'decodability_mean', 'decodability_std')
            allERPdata(count,:) = decodability_mean; allERPdata_std = decodability_std;
            allERPdata_sem(count,:) = allERPdata_std/nReps;
            
            loadnameNonPara = ['Results_Decode_FvR_NonPara_GC_a=' num2str(this_a) '_c=' num2str(this_c) '_d=' num2str(this_d) '.mat']
            load([filePathCoords loadnameNonPara])
            results_NP = results; clear results; clear data;
            
            
            %Decodability results (gather all)
            allPFC_logpwr(ff,vv) = results_NP.multiLogPwr_PFC_mean;
            allPFC_logpwr_std(ff,vv) = results_NP.multiLogPwr_PFC_std;
            allPFC_logpwr_sem(ff,vv) = allPFC_logpwr_std(ff,vv)/sqrt(length(results_NP.multiLogPwr_PFC));
            alllogPwr(ff,vv) = results_NP.multiLogPwr_FFA_mean;
            alllogPwr_std(ff,vv) = results_NP.multiLogPwr_FFA_std;
            alllogPwr_sem(ff,vv) = alllogPwr_std(ff,vv)/sqrt(length(results_NP.multiLogPwr_FFA));
            allCoh(ff,vv) = results_NP.multiCoh_mean;
            allCoh_sem(ff,vv) = results_NP.multiCoh_std;
            allCoh_sem(ff,vv) = allCoh_sem(ff,vv)/sqrt(length(results_NP.multiCoh));
            allGC_toPFC(ff,vv) = results_NP.multiGC_FFA_to_PFC_mean;
            allGC_toPFC_std(ff,vv) = results_NP.multiGC_FFA_to_PFC_std;
            allGC_toPFC_sem(ff,vv) = allGC_toPFC_std(ff,vv)/sqrt(length(results_NP.multiGC_FFA_to_PFC));
            allGC_fromPFC(ff,vv) = results_NP.multiGC_PFC_to_FFA_mean;
            allGC_fromPFC_std(ff,vv) = results_NP.multiGC_PFC_to_FFA_std;
            allGC_fromPFC_sem(ff,vv) = allGC_fromPFC_std(ff,vv)/sqrt(length(results_NP.multiGC_PFC_to_FFA));
        end
    end
end
%
%% Gather all data to plot
%GATHER ALL POWER DECODING RESULTS
allPFClogPwr = allPFC_logpwr(:)'; %PFC log power
allPFC_sem_logPwr = allPFC_logpwr_sem(:)';
allOtherlogPwr = alllogPwr(:)'; %Other node log power
allOther_sem_logPwr = alllogPwr_sem(:)';
alllogPwrs = [allPFClogPwr; allOtherlogPwr]; %Combined power
all_sem_logPwrs = [allPFC_sem_logPwr; allOther_sem_logPwr];

% Gather ALL ERP DECODING RESULTS
allERP = allERPdata'; % ERP results
allERP_sem = allERPdata_sem';

% Gather ALL Coherence DECODING RESULTS
allCoh_pts = allCoh(:)'; % Coherence
allCoh_sem_pts = allCoh_sem(:)';

%%Gather GC to PFC results
allGC_toPFC_results = allGC_toPFC(:)'; %GC to PFC
allGC_toPFC_sem_results = allGC_toPFC_sem(:)';

%%Gather GC from PFC results
allGC_fromPFC_results = allGC_fromPFC(:)'; %GC from PFC
allGC_fromPFC_sem_results = allGC_fromPFC_sem(:)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT ERP and log POWER (in one figure) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gather the data to plot
plot_ERP_dataALL = allERP;
plot_sem_ERP_dataALL = allERP_sem;
plot_pwr_ALL = alllogPwrs;
plot_sem_pwr_ALL = all_sem_logPwrs;

figure(2)
errorbar(plot_ERP_dataALL(1,:)',plot_sem_ERP_dataALL(1,:)','b-','LineWidth', 8); hold on;
errorbar(plot_ERP_dataALL(2,:)',plot_sem_ERP_dataALL(2,:)','k-','LineWidth', 8); hold on;
errorbar(plot_pwr_ALL(1,:)',plot_sem_pwr_ALL(1,:)','b-', 'LineWidth', 3); hold on;
errorbar(plot_pwr_ALL(2,:)',plot_sem_pwr_ALL(2,:)', 'k-', 'LineWidth', 3); hold on;
ylabel('Classification Accuracy (%)','FontSize',30);
set(gca, 'XTick', 1:length(plot_ERP_dataALL), 'XTickLabel', labels, 'FontSize',20);
xlim([0 length(plot_ERP_dataALL)+1]); ylim([45 100]);
title(['Decodability Waveforms and log Power - Face v Random Trials'], 'FontSize', 30);
plot(0:length(plot_ERP_dataALL)+1, ones(1,length(plot_ERP_dataALL)+2)*50, '-r', 'LineWidth', 2); hold on;
legend({'ERP PFC','ERP FFA','log Power PFC','log Power FFA'})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT COHERENCE DECODING RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gather the data to plot
allCoh_plot = allCoh_pts;
allCoh_sem = allCoh_sem_pts;

figure(22)
errorbar(allCoh_plot(1,:),allCoh_sem(1,:),'b-','LineWidth', 6); hold on;
ylabel('Classification Accuracy (%)','FontSize',30);
set(gca, 'XTick', 1:length(allCoh_plot), 'XTickLabel', labels, 'FontSize',20);
xlim([0 length(allCoh_plot)+1]); ylim([45 100]);
title(['Decodability Coherence - Face v Random Trials'], 'FontSize', 30);
plot(0:length(allCoh_plot)+1, ones(1,length(allCoh_plot)+2)*50, '-r', 'LineWidth', 2); hold on;
legend({'PFC and FFA'}, 'FontSize', 20)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT Granger Causality DECODING RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gather the data to plot
plot_GC_to_PFC_ALL = allGC_toPFC_results;
plot_sem_GC_to_PFC_ALL = allGC_toPFC_sem_results;
plot_GC_from_PFC_ALL = allGC_fromPFC_results;
plot_sem_GC_from_PFC_ALL = allGC_fromPFC_sem_results;

figure(212)
errorbar(plot_GC_from_PFC_ALL(1,:),plot_sem_GC_from_PFC_ALL(1,:),'LineWidth', 5); hold on;
errorbar(plot_GC_to_PFC_ALL(1,:),plot_sem_GC_to_PFC_ALL(1,:),'LineWidth', 5); hold on;
ylabel('Classification Accuracy (%)','FontSize',30);
set(gca, 'XTick', 1:length(plot_GC_from_PFC_ALL), 'XTickLabel', labels, 'FontSize',20);
xlim([0 length(plot_GC_from_PFC_ALL)+1]); ylim([45 100]);
title(['Decodability Granger Causality - Face v Random Trials'], 'FontSize', 30);
plot(0:length(plot_GC_from_PFC_ALL)+1, ones(1,length(plot_GC_from_PFC_ALL)+2)*50, '-r', 'LineWidth', 2); hold on;
legend({'PFC to FFA', 'FFA to PFC'})
