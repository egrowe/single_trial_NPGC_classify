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

filePathCoords = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/3_Decoding/'];

%A matrix connectivity
poss_a = -0.5: 0.05: 0.5; %weight for 'a' (= x self conn X)
poss_c = -0.5: 0.05: 0.5; %weight for 'c' (= x self conn X in Y data)
poss_d = -0.5: 0.05: 0.5;  %weight for 'd' (= y self conn Y)
nReps = 10; %repeat decoding scheme this many times

count = 0; countA = 1; countB = 1; countC = 1; %setup indexes

for hh = 1:length(poss_a)
    
    for  ff = 1:length(poss_c)

        for vv = 1:length(poss_d)
            count = count + 1;

            %Model parameters
            this_a = poss_a(gg); % use this weight for 'a' (= x self conn X)
            this_c = poss_c(ff);  % use this weight for 'c' (= x self conn X in Y data)
            this_d = poss_d(vv); % use this weight for 'd' (= y self conn Y)

            %Load this datafile
            loadnameNonPara = ['Results_Decode_FvR_NonPara_GC_a=' num2str(this_a) '_c=' num2str(this_c) '_d=' num2str(this_d) '.mat'];
            load([filePathCoords loadnameNonPara])
            
            data_NP = data; clear results; clear data; %ensure old data is overwritten
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% ANALYSIS 1a: COHERENCE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cohAll_faces = [];  cohAll_random = []; %create empty variables
            %% Load Face Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            freqWindow = linspace(0,40,20)
            
            cohAll_faces =  data_NP.coh_faces'% extract only value that matters
            
            %Find the mean coherence value across all trials (per frequency)
            mCoh_faces = squeeze(mean(cohAll_faces,1)); %mean all coh vals
            stdCoh_faces = squeeze(std(cohAll_faces,1)); %std all coh vals
            semCoh_faces = stdCoh_faces/(sqrt(length(cohAll_faces))); % SEM coherence between ROIs
            
            %Setup plot indexes
            usePlotIdxs = [1:18];
            
            %PLOT COHERENCE VALUES: Faces
            figure(1)
            subplot(2,4,countA)
            errorbar(mCoh_faces,semCoh_faces, 'LineWidth', 3); hold on; ylim([0.35 0.9]);
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {ceil(freqWindow(1)),ceil(freqWindow(10)),ceil(freqWindow(20))}); xlim([1 length(mCoh_faces)]);
            % set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]); %end %set(gca, 'FontSize', 1); end
            %title(['COHERENCE: ' num2str(useMainCoord) ' and ' num2str(useSupCoord)], 'FontSize', 15); %legend({'', 'TEST', '', 'TRAIN'});
            hold on;
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Load Random Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cohAll_random = data_NP.coh_random'; %cohAll_random(:,:,2); % extract only value that matters
            
            %Find the mean coherence value across all trials (per frequency)
            mCoh_random = squeeze(mean(cohAll_random,1)); %mean all coh vals
            stdCoh_random = squeeze(std(cohAll_random,1)); %std all coh vals
            semCoh_random = stdCoh_random/(sqrt(length(cohAll_random))); % SEM coherence between ROIs
            
            figure(1)
            subplot(2,4,countA)
            errorbar(mCoh_random,semCoh_random, 'LineWidth', 3); hold on;
            countA = countA + 1;
            legend({'FACES', 'RANDOM'});
            
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
            
            
            
            GO_FROM_PFC_faces = mGC_faces_FROM_PFC;
            GC_sem_FROM_PFC_faces = semGC_faces_FROM_PFC;
            
            figure(333)
            subplot(2,4,countC)
            errorbar(GC_TO_PFC_faces,GC_sem_TO_PFC_faces, 'LineWidth', 3); hold on; ylim([0.1 0.9]);
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {ceil(freqWindow(1)),ceil(freqWindow(10)),ceil(freqWindow(20))}); xlim([1 length(mCoh_faces)]);
            % set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]); %set(gca, 'FontSize', 1); end
            % title(['GC to ' num2str(useMainCoord) ' from:' num2str(useSupCoord)], 'FontSize', 15); %legend({'', 'TEST', '', 'TRAIN'});
            hold on;
            
            figure(444)
            %thisPlot = usePlotIdxs(countA);
            subplot(2,4,countC)
            errorbar(GO_FROM_PFC_faces,GC_sem_FROM_PFC_faces, 'LineWidth', 3); hold on; ylim([0.1 0.9]);
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {ceil(freqWindow(1)),ceil(freqWindow(10)),ceil(freqWindow(20))}); xlim([1 length(mCoh_faces)]);% else
            %   set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]); end %set(gca, 'FontSize', 1); end
            %title(['GC from ' num2str(useMainCoord) ' to: ' num2str(useSupCoord)], 'FontSize', 15); %legend({'', 'TEST', '', 'TRAIN'});
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
            
            GO_FROM_PFC_random = mGC_random_FROM_PFC;%(:,2)
            GC_sem_FROM_PFC_random = semGC_random_FROM_PFC;%(:,2)
            
            figure(333)
            subplot(2,4,countC)
            errorbar(GC_TO_PFC_random,GC_sem_TO_PFC_random, 'LineWidth', 3); hold on;
            
            figure(444)
            subplot(2,4,countC)
            errorbar(GO_FROM_PFC_random,GC_sem_FROM_PFC_random, 'LineWidth', 3); hold on;
            
            countC = countC + 1;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% ANALYSIS 2a: POWER RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %% Load Face Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            load(loadNameFaces); %load face data
            
            for ss = 1:size(allPowerInfo,1)
                pwrData = allPowerInfo(ss).powspctrm(:,1:length(freqWindow))
                pwrAll_faces(ss,:,:) = pwrData'; % power (form trials x freq x pwrPerCond)
            end
            
            %Extract log power for FFA
            mlogPwr_FFA_faces = squeeze(mean(data_NP.Logpwr_FFA_faces,2)'); %MEAN power across trials
            stdlogPwr_FFA_faces = squeeze(std(data_NP.Logpwr_FFA_faces')); % STDEV power across trials
            semlogPwr_FFA_faces = stdlogPwr_FFA_faces/(sqrt(length(data_NP.Logpwr_FFA_faces))); %SEM power across trials
            
            %Extract log power for PFC
            mlogPwr_PFC_faces = squeeze(mean(data_NP.Logpwr_PFC_faces,2)'); %MEAN power across trials
            stdlogPwr_PFC_faces = squeeze(std(data_NP.Logpwr_PFC_faces')); % STDEV power across trials
            semlogPwr_PFC_faces = stdlogPwr_PFC_faces/(sqrt(length(data_NP.Logpwr_PFC_faces))); %SEM power across trials
            
            
            
            figure(1000)
            thisPlot = usePlotIdxs(countB);
            subplot(2,4,countB)
            %Plot for FFA
            errorbar(mlogPwr_FFA_faces,semlogPwr_FFA_faces, 'LineWidth', 2); hold on;
            if countB == 1
                ylim([-2.1 0.2]); xlim([1 length(mlogPwr_FFA_faces)]); % ylim([-3 0]); %ylabel('Coherence','FontSize',12);
                set(gca, 'XTick', [1,10,21], 'XTickLabel', {ceil(freqWindow(1)),ceil(freqWindow(10)),ceil(freqWindow(20))});
                set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15); xlim([1 length(mCoh_faces)]);% else
            else set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]); xlim([1 length(mCoh_faces)]); end
            
            %Plot for PFC
            subplot(2,4,countB+4)
            
            errorbar(mlogPwr_PFC_faces,semlogPwr_PFC_faces, 'LineWidth', 2); hold on;
            if countB == 1
                ylim([-2.1 0.2]); xlim([1 length(mlogPwr_FFA_faces)]); % ylim([-3 0]); %ylabel('Coherence','FontSize',12);
                set(gca, 'XTick', [1,10,21], 'XTickLabel', {ceil(freqWindow(1)),ceil(freqWindow(10)),ceil(freqWindow(20))});
                set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15); xlim([1 length(mCoh_faces)]);% else
            else set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]); xlim([1 length(mCoh_faces)]); end
            
            
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Load Random Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Extract log power for FFA
            mlogPwr_FFA_random = squeeze(mean(data_NP.Logpwr_FFA_random,2)'); %MEAN power across trials
            stdlogPwr_FFA_random = squeeze(std(data_NP.Logpwr_FFA_random')); % STDEV power across trials
            semlogPwr_FFA_random = stdlogPwr_FFA_random/(sqrt(length(data_NP.Logpwr_FFA_random))); %SEM power across trials
            
            %Extract log power for PFC
            mlogPwr_PFC_random = squeeze(mean(data_NP.Logpwr_PFC_random,2)'); %MEAN power across trials
            stdlogPwr_PFC_random = squeeze(std(data_NP.Logpwr_PFC_random')); % STDEV power across trials
            semlogPwr_PFC_random = stdlogPwr_PFC_random/(sqrt(length(data_NP.Logpwr_PFC_random))); %SEM power across trials
            
            figure(1000)
            subplot(2,4,countB)
            %Plot for FFA
            errorbar(mlogPwr_FFA_random,semlogPwr_FFA_random, 'LineWidth', 2); hold on;
            
            %Plot for PFC
            subplot(2,4,countB+4)
            errorbar(mlogPwr_PFC_random,semlogPwr_PFC_random, 'LineWidth', 2); hold on;
            countB = countB + 1;
            
            
        end
    end
end