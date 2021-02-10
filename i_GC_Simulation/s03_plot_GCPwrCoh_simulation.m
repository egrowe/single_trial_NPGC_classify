%% Visualise NP coherence, power and Granger Causality results - SIMULATION

% This script enable visualisation of the results from non-parametric
% Granger Causality analysis including coherence, power and GC
% coefficients.

% This script requires the libSVM toolbox is added to the path

% Each stage split into VISUALISATION raw data, DECODING data and
% visualisation of the decoding results

% Extract data for visualisation and decoding for either:
% (1): Coherence
% (2): Log Power
% (3): Granger causality


%% PLOT ALL SIMULATION RESULTS
clear all; close all;
addpath('/Users/egrow1/Desktop/Shafto_GitHub_Scripts/i_GC_Simulation/mseb');
filePathCoords = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/i_GC_Simulation/'];

%A matrix connectivity
poss_a = 0.5% [-0.5: 0.5: 0.5]; %weight for 'a' (= x self conn X)
poss_c = 0.5%[-0.5: 0.5: 0.5]; %weight for 'c' (= x self conn X in Y data)
poss_d = 0.5%[-0.5: 0.5: 0.5];  %weight for 'd' (= y self conn Y)

count = 0; countA = 1; countB = 1; countC = 1; %setup indexes

for gg = 1:length(poss_a)
    
    for  ff = 1:length(poss_c)

        for vv = 1:length(poss_d)
            count = count + 1;

            %Model parameters
            this_a = poss_a(gg); % use this weight for 'a' (= x self conn X)
            this_c = poss_c(ff);  % use this weight for 'c' (= x self conn X in Y data)
            this_d = poss_d(vv); % use this weight for 'd' (= y self conn Y)

            %Load this datafile
   %         loadnameNonPara = ['Results_Decode_FvR_NonPara_GC_XX=' num2str(this_a) '_XY=' num2str(this_c) '_YY=' num2str(this_d) '.mat'];
             loadnameNonPara = (['Results_Decode_55Hz_FvR_NonPara_GC_XX=' num2str(this_a) '_XY=' num2str(this_c) '_YY=' num2str(this_d) '.mat']);
             load([filePathCoords loadnameNonPara])

            data_NP = data; results_NP = results; clear results; clear data; %ensure old data is overwritten
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% ANALYSIS 1a: COHERENCE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cohAll_faces = [];  cohAll_random = []; %create empty variables
            %% Load Face Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            freqWindow = linspace(0,53.7849,28)
            
            cohAll_faces =  data_NP.coh_faces'% extract only value that matters
            
            %Find the mean coherence value across all trials (per frequency)
            mCoh_faces = squeeze(mean(cohAll_faces,1)); %mean all coh vals
            stdCoh_faces = squeeze(std(cohAll_faces,1)); %std all coh vals
            semCoh_faces = stdCoh_faces/(sqrt(length(cohAll_faces))); % SEM coherence between ROIs
            
            %PLOT COHERENCE VALUES: Faces
            figure(1)
            subplot(3,9,countA)
            errorbar(mCoh_faces,semCoh_faces, 'LineWidth', 3); hold on; ylim([0.35 0.9]);
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),ceil(length(freqWindow))], 'XTickLabel', {ceil(freqWindow(1)),freqWindow(ceil(length(freqWindow)/2)),freqWindow(ceil(length(freqWindow)))}); xlim([1 length(mCoh_faces)]);
            title(['XX: ' num2str(this_a) ' XY:' num2str(this_c) ' YY:'  num2str(this_d) ], 'FontSize', 10); %legend({'', 'TEST', '', 'TRAIN'});
            hold on;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Load Random Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cohAll_random = data_NP.coh_random'; %cohAll_random(:,:,2); % extract only value that matters
            
            %Find the mean coherence value across all trials (per frequency)
            mCoh_random = squeeze(mean(cohAll_random,1)); %mean all coh vals
            stdCoh_random = squeeze(std(cohAll_random,1)); %std all coh vals
            semCoh_random = stdCoh_random/(sqrt(length(cohAll_random))); % SEM coherence between ROIs
            
            figure(1)
            subplot(3,9,countA)
            errorbar(mCoh_random,semCoh_random, 'LineWidth', 3); hold on;
            %legend({'FACES', 'RANDOM'});
            
            %Coh data for face and random
            all_m_Coh = [mCoh_faces; mCoh_random]
            all_sem_Coh = [semCoh_faces; semCoh_random]

          %  figure(11111)
          %  subplot(3,9,countA)
                     figure(88)
           subplot(1,5,3)
            mseb((repmat(1:length(all_m_Coh),2,1)),all_m_Coh,all_sem_Coh,[],1)
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),ceil(length(freqWindow))], 'XTickLabel', {ceil(freqWindow(1)),freqWindow(ceil(length(freqWindow)/2)),freqWindow(ceil(length(freqWindow)))}); xlim([1 length(mCoh_faces)]);
            title(['XX: ' num2str(this_a) ' XY:' num2str(this_c) ' YY:'  num2str(this_d) ], 'FontSize', 10); %legend({'', 'TEST', '', 'TRAIN'});
            
            countA = countA + 1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% ANALYSIS 3a: GRANGER CAUSALITY RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %
            GCall_faces_TO_PFC = data_NP.GC_FFA_to_PFC_faces';
            
            %Gather mean and std of GC
            mGC_faces_TO_PFC = squeeze(mean(GCall_faces_TO_PFC,1)); %MEAN GC over trials (at each frequency band)
            stdGC_faces_TO_PFC = squeeze(std(GCall_faces_TO_PFC)); %STD GC over trials
            semGC_faces_TO_PFC = stdGC_faces_TO_PFC/(sqrt(length(GCall_faces_TO_PFC))); %SEM GC over trials
            
            
            GC_TO_PFC_faces = mGC_faces_TO_PFC;
            GC_sem_TO_PFC_faces = semGC_faces_TO_PFC;
            
            
            GC_all_faces_FROM_PFC = data_NP.GC_PFC_to_FFA_faces';
            
            %Gather mean and std of GC
            mGC_faces_FROM_PFC = squeeze(mean(GC_all_faces_FROM_PFC,1)); %MEAN GC over trials (at each frequency band)
            stdGC_faces_FROM_PFC = squeeze(std(GC_all_faces_FROM_PFC)); %STD GC over trials
            semGC_faces_FROM_PFC = stdGC_faces_FROM_PFC/(sqrt(length(GC_all_faces_FROM_PFC))); %SEM GC over trials

            
            GC_FROM_PFC_faces = mGC_faces_FROM_PFC;
            GC_sem_FROM_PFC_faces = semGC_faces_FROM_PFC;
            
            figure(333)
            subplot(3,9,countC)
            errorbar(GC_TO_PFC_faces,GC_sem_TO_PFC_faces, 'LineWidth', 3); hold on; ylim([0.1 1]);
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),ceil(length(freqWindow))], 'XTickLabel', {ceil(freqWindow(1)),freqWindow(ceil(length(freqWindow)/2)),freqWindow(ceil(length(freqWindow)))}); xlim([1 length(mCoh_faces)]);
            title(['XX: ' num2str(this_a) ' XY:' num2str(this_c) ' YY:'  num2str(this_d) ], 'FontSize', 10); %legend({'', 'TEST', '', 'TRAIN'});
            
            figure(444)
            subplot(3,9,countC)
            errorbar(GC_FROM_PFC_faces,GC_sem_FROM_PFC_faces, 'LineWidth', 3); hold on; ylim([0.1 1]);
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),ceil(length(freqWindow))], 'XTickLabel', {ceil(freqWindow(1)),freqWindow(ceil(length(freqWindow)/2)),freqWindow(ceil(length(freqWindow)))}); xlim([1 length(mCoh_faces)]);
            title(['XX: ' num2str(this_a) ' XY:' num2str(this_c) ' YY:'  num2str(this_d) ], 'FontSize', 10); %legend({'', 'TEST', '', 'TRAIN'});
            hold on;
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Load Random Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            GCall_random_TO_PFC = data_NP.GC_FFA_to_PFC_random';
            
            %Gather mean and std of GC
            mGC_random_TO_PFC = squeeze(mean(GCall_random_TO_PFC,1)); %MEAN GC over trials (at each frequency band)
            stdGC_random_TO_PFC = squeeze(std(GCall_random_TO_PFC)); %STD GC over trials
            semGC_random_TO_PFC = stdGC_random_TO_PFC/(sqrt(length(GCall_random_TO_PFC))); %SEM GC over trials
            
            GC_TO_PFC_random = mGC_random_TO_PFC;
            GC_sem_TO_PFC_random = semGC_random_TO_PFC;
            
            GC_all_random_FROM_PFC = data_NP.GC_PFC_to_FFA_random';
            
            %Gather mean and std of GC
            mGC_random_FROM_PFC = squeeze(mean(GC_all_random_FROM_PFC,1)); %MEAN GC over trials (at each frequency band)
            stdGC_random_FROM_PFC = squeeze(std(GC_all_random_FROM_PFC)); %STD GC over trials
            semGC_random_FROM_PFC = stdGC_random_FROM_PFC/(sqrt(length(GC_all_random_FROM_PFC))); %SEM GC over trials
            
            GC_FROM_PFC_random = mGC_random_FROM_PFC;%(:,2)
            GC_sem_FROM_PFC_random = semGC_random_FROM_PFC;%(:,2)
            
            figure(333)
            subplot(3,9,countC)
            errorbar(GC_TO_PFC_random,GC_sem_TO_PFC_random, 'LineWidth', 3); hold on;
            
            figure(444)
            subplot(3,9,countC)
            errorbar(GC_FROM_PFC_random,GC_sem_FROM_PFC_random, 'LineWidth', 3); hold on;

            %Gather both data
            all_GC_to_PFC = [GC_TO_PFC_faces; GC_TO_PFC_random];
            all_GC_to_PFC_sem = [GC_sem_TO_PFC_faces; GC_sem_TO_PFC_random];
            
            all_GC_from_PFC = [GC_FROM_PFC_faces; GC_FROM_PFC_random];
            all_GC_from_PFC_sem = [GC_sem_FROM_PFC_faces; GC_sem_FROM_PFC_random];
            
            
        %    figure(31111)
        %    subplot(3,9,countC)
                   figure(88)
           subplot(1,5,4)
            mseb((repmat(1:length(all_GC_to_PFC),2,1)),all_GC_to_PFC,all_GC_to_PFC_sem,[],1)
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),ceil(length(freqWindow))], 'XTickLabel', {ceil(freqWindow(1)),...
                freqWindow(ceil(length(freqWindow)/2)),freqWindow(ceil(length(freqWindow)))}); xlim([1 length(mCoh_faces)]);
            title(['XX: ' num2str(this_a) ' XY:' num2str(this_c) ' YY:'  num2str(this_d) ], 'FontSize', 10); 
            
          %  figure(41111)
         %   subplot(3,9,countC)
                    figure(88)
           subplot(1,5,5)
            mseb((repmat(1:length(all_GC_from_PFC),2,1)),all_GC_from_PFC,all_GC_from_PFC_sem,[],1)
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),ceil(length(freqWindow))], 'XTickLabel', {ceil(freqWindow(1)),...
                freqWindow(ceil(length(freqWindow)/2)),freqWindow(ceil(length(freqWindow)))}); xlim([1 length(mCoh_faces)]);
            title(['XX: ' num2str(this_a) ' XY:' num2str(this_c) ' YY:'  num2str(this_d) ], 'FontSize', 10); 

            countC = countC + 1;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% ANALYSIS 2a: POWER RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %% Load Face Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            
            %Extract log power for FFA
            Pwr_FFA_faces = data_NP.Logpwr_FFA_faces; %MEAN power across trials
            log_Pwr_FFA_faces = log10(Pwr_FFA_faces+10);
            mlogPwr_FFA_faces = mean(log_Pwr_FFA_faces,2)';
            stdlogPwr_FFA_faces = squeeze(std(log_Pwr_FFA_faces')); % STDEV power across trials
            semlogPwr_FFA_faces = stdlogPwr_FFA_faces/(sqrt(length(log_Pwr_FFA_faces))); %SEM power across trials
            
            %Extract log power for PFC
            Pwr_PFC_faces = data_NP.Logpwr_PFC_faces;
            log_Pwr_PFC_faces = log10(Pwr_PFC_faces+10);
            mlogPwr_PFC_faces = squeeze(mean(log_Pwr_PFC_faces,2)'); %MEAN power across trials
            stdlogPwr_PFC_faces = squeeze(std(log_Pwr_PFC_faces')); % STDEV power across trials
            semlogPwr_PFC_faces = stdlogPwr_PFC_faces/(sqrt(length(log_Pwr_PFC_faces))); %SEM power across trials
            
            figure(1000)
            subplot(3,9,countB)
            %Plot for FFA
            errorbar(mlogPwr_FFA_faces,semlogPwr_FFA_faces, 'LineWidth', 2); hold on;
            if countB == 1
            set(gca, 'XTick', [1,round(length(freqWindow)/2),ceil(length(freqWindow))], ...
                'XTickLabel', {ceil(freqWindow(1)),freqWindow(ceil(length(freqWindow)/2)),freqWindow(ceil(length(freqWindow)))}); xlim([1 length(mCoh_faces)]);
                set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15); xlim([1 length(mCoh_faces)]); ylim([0.8 0.95]);% else
            else set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]); xlim([1 length(mCoh_faces)]); ylim([0.8 0.95]); end
            title(['XX: ' num2str(this_a) ' XY:' num2str(this_c) ' YY:'  num2str(this_d) ], 'FontSize', 10); hold on;

            %Plot for PFC
            figure(1001)
            subplot(3,9,countB)
            errorbar(mlogPwr_PFC_faces,semlogPwr_PFC_faces, 'LineWidth', 2); hold on;
            if countB == 1
            set(gca, 'XTick', [1,round(length(freqWindow)/2),ceil(length(freqWindow))], 'XTickLabel', {ceil(freqWindow(1)),freqWindow(ceil(length(freqWindow)/2)),freqWindow(ceil(length(freqWindow)))}); xlim([1 length(mCoh_faces)]);
                set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15); xlim([1 length(mCoh_faces)]);
            else set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]); xlim([1 length(mCoh_faces)]); ylim([0.8 0.95]); end
            ylim([0.8 0.95]); xlim([1 length(mlogPwr_FFA_faces)]);
            title(['XX: ' num2str(this_a) ' XY:' num2str(this_c) ' YY:'  num2str(this_d) ], 'FontSize', 10); hold on;
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Load Random Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Extract log power for FFA
            Pwr_FFA_random = data_NP.Logpwr_FFA_random;
            log_Pwr_FFA_random = log10(Pwr_FFA_random+10)
            mlogPwr_FFA_random = mean(log_Pwr_FFA_random,2)'; %MEAN power across trials
            stdlogPwr_FFA_random = squeeze(std(log_Pwr_FFA_random')); % STDEV power across trials
            semlogPwr_FFA_random = stdlogPwr_FFA_random/(sqrt(length(log_Pwr_FFA_random))); %SEM power across trials
            
            
            %Extract log power for PFC
            Pwr_PFC_random = data_NP.Logpwr_PFC_random; %MEAN power across trials
            log_Pwr_PFC_random = log10(Pwr_PFC_random+10);
            mlogPwr_PFC_random = mean(log_Pwr_PFC_random,2)'
            stdlogPwr_PFC_random = squeeze(std(log_Pwr_PFC_random')); % STDEV power across trials
            semlogPwr_PFC_random = stdlogPwr_PFC_random/(sqrt(length(log_Pwr_PFC_random))); %SEM power across trials
            
            figure(1000)
            subplot(3,9,countB)
            %Plot for FFA
            errorbar(mlogPwr_FFA_random,semlogPwr_FFA_random, 'LineWidth', 2); hold off;
            
            %Plot for PFC
            figure(1001)
            subplot(3,9,countB)
            errorbar(mlogPwr_PFC_random,semlogPwr_PFC_random, 'LineWidth', 2); hold off;
            
            %All power data
            all_Logpwr_PFC = [mlogPwr_PFC_faces; mlogPwr_PFC_random]
            all_Logpwr_PFC_sem = [semlogPwr_PFC_faces; semlogPwr_PFC_random]
            all_Logpwr_FFA = [mlogPwr_FFA_faces; mlogPwr_FFA_random]
            all_Logpwr_FFA_sem = [semlogPwr_FFA_faces; semlogPwr_FFA_random]
            
           % figure(10022)
           % subplot(3,9,countB)
           figure(88)
           subplot(1,5,1)
         
            %Plot for FFA
            mseb((repmat(1:length(all_Logpwr_FFA),2,1)),all_Logpwr_FFA,all_Logpwr_FFA_sem,[],1)
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),ceil(length(freqWindow))], 'XTickLabel', {ceil(freqWindow(1)),...
                freqWindow(ceil(length(freqWindow)/2)),freqWindow(ceil(length(freqWindow)))}); xlim([1 length(mCoh_faces)]);
            title(['XX: ' num2str(this_a) ' XY:' num2str(this_c) ' YY:'  num2str(this_d) ], 'FontSize', 10); 

           % figure(10033)
          %  subplot(3,9,countB)
           figure(88)
           subplot(1,5,2)
            %Plot for FFA
            mseb((repmat(1:length(all_Logpwr_PFC),2,1)),all_Logpwr_PFC,all_Logpwr_PFC_sem,[],1)
          set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),ceil(length(freqWindow))], 'XTickLabel', {ceil(freqWindow(1)),...
                freqWindow(ceil(length(freqWindow)/2)),freqWindow(ceil(length(freqWindow)))}); xlim([1 length(mCoh_faces)]);
            title(['XX: ' num2str(this_a) ' XY:' num2str(this_c) ' YY:'  num2str(this_d) ], 'FontSize', 10); 

            countB = countB + 1;
            
            
        end
    end
end

%% DECODE AND PLOT