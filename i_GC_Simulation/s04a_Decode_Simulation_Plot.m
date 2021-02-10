%% Extract FACES coherence, power and GC for single-trial estimates
close all; clear all;
% Add the decoding libSVM folder to the path
addpath(genpath('/Users/egrow1/Desktop/NEW_Redo_OriginalChans_SHAFTO/decode_libSVM'))
addpath(genpath('/Users/egrow1/Desktop/Shafto_GitHub_Scripts/3_Decoding'))

%Set the filepath to load data from
filepath = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/i_GC_Simulation/'];

%A matrix connectivity
poss_a = 0.5 %[0.5]; %weight for 'a' (= x self conn X)
poss_c = 0.5% [-0.5, 0, 0.5]; %weight for 'c' (= x self conn X in Y data)
poss_d = 0.5%[-0.5, 0, 0.5];  %weight for 'd' (= y self conn Y)
maxFreq = 55;

%Decoding parameters
zscore = 1; % set zscore on (decoding)
nReps = 10; %number of times to repeat decoding scheme
count = 0;  hh = 1; %setup indexes

% FIRST OF ALL, DECODE USING SOURCE WAVEFORMS
%Load name for files
loadname = ([filepath 'All_Data_Synthetic_SourceERP_multi_Xs=-0.5_XsY=0.5_Ys=0.5.mat'])
load(loadname)
% DECODE USING THE WAVEFORMS
costRange = -5;
[decodability_sWaves_nodeY, cv_acc, w] = decode_multiFeature_libSVM(YF_all', YR_all', nReps, costRange, zscore)
[decodability_sWaves_nodeX, cv_acc, w] = decode_multiFeature_libSVM(XF_all', XR_all', nReps, costRange, zscore)


for hh = 1:length(poss_a)
    
    for  ff = 1:length(poss_c)
        
        for vv = 1:length(poss_d)
            count = count + 1;
            
            %Model parameters
            this_a = poss_a(hh); % use this weight for 'a' (= x self conn X)
            this_c = poss_c(ff);  % use this weight for 'c' (= x self conn X in Y data)
            this_d = poss_d(vv); % use this weight for 'd' (= y self conn Y)
            
            load(['allfaces_CohPower_GC_SYNTHETIC_XX=' num2str(this_a) '_XY=' num2str(this_c) '_YY=' num2str(this_d) '.mat'])
            
            %Save name
            savename = ([filepath 'Results_Decode_55Hz_FvR_NonPara_GC_XX=' num2str(this_a) '_XY=' num2str(this_c) '_YY=' num2str(this_d) '.mat']);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% COHERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %First look at coherence per trial
            freqs =  allGranger(1).freq
            freq_roi = 0 <= freqs & freqs <= maxFreq
            freqWindow = freqs(freq_roi)
            freqIdxs = 1:length(freqWindow)
            
            for trialNo = 1:size(allCoherenceInfo,1)
                for thisFreq = 1:length(freqWindow)
                    cohData = allCoherenceInfo(trialNo).cohspctrm(:,:,thisFreq)
                    cohAll(trialNo,thisFreq,:) = cohData(:)'; % coherence (form trials x freq x cohVals)
                end
            end
            
            %Compute mean and std COHERENCE
            mean_Coh = squeeze(mean(cohAll,1))
            std_Coh = squeeze(std(cohAll))
            mean_Coh_b = mean_Coh(:,2)
            std_Coh_b = std_Coh(:,2)
            sem_Coh_b = std_Coh_b/(sqrt(length(std_Coh_b)))
            
            %PLOT
            figure(1)
            errorbar(mean_Coh_b,sem_Coh_b, 'LineWidth', 4); hold on;
            set(gca, 'FontSize', 30); xlabel('Frequency (Hz)','FontSize',30);
            ylabel('Coherence (mean across trials)','FontSize',20);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {freqWindow(1),freqWindow(10),freqWindow(20)}); xlim([1 length(mean_Coh_b)]);
            title('COHERENCE: Left PFC and left FFA', 'FontSize', 20); %legend({'', 'TEST', '', 'TRAIN'});
            hold on;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% POWER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for ss = 1:size(allPowerInfo,1)
                pwrData = allPowerInfo(ss).powspctrm(:,1:length(freqWindow))
                pwr_All(ss,:,:) = pwrData'; % power (form trials x freq x pwrPerCond)
            end
            
            mean_pwr = squeeze(mean(pwr_All,1))
            std_pwr = squeeze(std(pwr_All))
            sem_pwr = std_pwr/(sqrt(length(std_pwr)))
            
            %Log power
            logPwrAll = log10(pwr_All);
            mean_log_pwr = squeeze(mean(logPwrAll,1))
            std_log_pwr = squeeze(std(logPwrAll))
            sem_log_GC = std_log_pwr/(sqrt(length(std_log_pwr)))
            
            
            figure(21)
            errorbar(mean_pwr(:,1),sem_pwr(:,1), 'LineWidth', 4); hold on;
            set(gca, 'FontSize', 30); xlabel('Frequency (Hz)','FontSize',30);
            ylabel('Power (mean across trials)','FontSize',30);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {freqWindow(1),freqWindow(10),freqWindow(20)});
            xlim([1 length(mean_Coh_b)]); title('POWER - Left  FFA', 'FontSize', 20); hold on;
            
            figure(22)
            errorbar(mean_pwr(:,2),sem_pwr(:,2), 'LineWidth', 4); hold on;
            set(gca, 'FontSize', 30); xlabel('Frequency (Hz)','FontSize',30);
            ylabel('Power (mean across trials)','FontSize',30);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {freqWindow(1),freqWindow(10),freqWindow(20)});
            xlim([1 length(mean_Coh_b)]); title('POWER - Left PFC', 'FontSize', 20); hold on;
            
            figure(23)
            imagesc(mean_pwr'); hold on;
            set(gca, 'FontSize', 30); xlabel('Frequency (Hz)','FontSize',30);
            ylabel('Power (mean across trials)','FontSize',30);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {freqWindow(1),freqWindow(10),freqWindow(20)}, 'YTick',[1,2],'YTickLabel',{'FFA','PFC'});
            xlim([1 length(mean_Coh_b)]); title('POWER - Left  FFA and PFC - FACE TRIALS', 'FontSize', 20); hold on;
            
            
            figure(24)
            errorbar(mean_log_pwr(:,1),sem_log_GC(:,1), 'LineWidth', 4); hold on;
            set(gca, 'FontSize', 30); xlabel('Frequency (Hz)','FontSize',30);
            ylabel('Log power (mean across trials)','FontSize',30);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {freqWindow(1),freqWindow(10),freqWindow(20)});
            xlim([1 length(mean_Coh_b)]); title('LOG POWER - Left  FFA', 'FontSize', 20); hold on;
            
            figure(25)
            errorbar(mean_log_pwr(:,2),sem_log_GC(:,2), 'LineWidth', 4); hold on;
            set(gca, 'FontSize', 30); xlabel('Frequency (Hz)','FontSize',30);
            ylabel('Log power (mean across trials)','FontSize',30);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {freqWindow(1),freqWindow(10),freqWindow(20)});
            xlim([1 length(mean_Coh_b)]); title('LOG POWER - Left PFC', 'FontSize', 20); hold on;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% GC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for trialNo = 1:size(allGranger,1)
                for thisFreq = 1:length(freqWindow)
                    GCData = allGranger(trialNo).grangerspctrm(:,:,thisFreq)
                    GCAll(trialNo,thisFreq,:) = GCData(:)'; % Granger (form trials x freq x GCVals)
                end
            end
            
            %Gather mean and std of GC
            mean_GC = squeeze(mean(GCAll,1))
            std_GC = squeeze(std(GCAll))
            sem_GC = std_GC/(sqrt(length(std_GC)))
            
            %Log GC
            useGC = GCAll(:,:,2:3)
            log_GCAll = log10(useGC)
            mean_log_GC = squeeze(mean(log_GCAll,1))
            std_log_GC = squeeze(std(log_GCAll))
            sem_log_GC = std_log_GC/(sqrt(length(std_log_GC)))
            
            figure(31)
            errorbar(mean_GC(:,2),sem_GC(:,2), 'LineWidth', 4); hold on;
            set(gca, 'FontSize', 30); xlabel('Frequency (Hz)','FontSize',30);
            ylabel('Granger Causality (mean across trials)','FontSize',20);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {freqWindow(1),freqWindow(10),freqWindow(20)}); xlim([1 length(mean_Coh_b)]);
            title('GC: Left PFC to Left FFA', 'FontSize', 20);
            
            figure(32)
            errorbar(mean_GC(:,3),sem_GC(:,3), 'LineWidth', 4); hold on;
            set(gca, 'FontSize', 30); xlabel('Frequency (Hz)','FontSize',30);
            ylabel('Granger Causality (mean across trials)','FontSize',20);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {freqWindow(1),freqWindow(10),freqWindow(20)}); xlim([1 length(mean_Coh_b)]);
            title('GC: Left FFA to Left PFC', 'FontSize', 20);
            
            figure(33)
            errorbar(mean_log_GC(:,1),sem_log_GC(:,1), 'LineWidth', 4); hold on;
            set(gca, 'FontSize', 30); xlabel('Frequency (Hz)','FontSize',30);
            ylabel('Log Granger Causality (mean across trials)','FontSize',20);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {freqWindow(1),freqWindow(10),freqWindow(20)}); xlim([1 length(mean_Coh_b)]);
            title('Log GC: Left PFC to Left FFA', 'FontSize', 20);
            
            figure(34)
            errorbar(mean_log_GC(:,2),sem_log_GC(:,2), 'LineWidth', 4); hold on;
            set(gca, 'FontSize', 30); xlabel('Frequency (Hz)','FontSize',30);
            ylabel('Log Granger Causality (mean across trials)','FontSize',20);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {freqWindow(1),freqWindow(10),freqWindow(20)}); xlim([1 length(mean_Coh_b)]);
            title('Log GC: Left FFA to Left PFC', 'FontSize', 20);
            
            plot_logGC_first_faces = log_GCAll(:,1,1);
            plot_GC_first_faces = GCAll(:,1,2);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% Get data for decoding
            
            % % DECODE COHERENCE
            coh_Decode_faces = cohAll(:,:,3)'; % want in form of features (coh at freq) x trials
            
            % DECODE POWER
            pwr_FFA_Decode_faces = pwr_All(:,:,1)';
            pwr_PFC_Decode_faces = pwr_All(:,:,2)';
            
            logPwr_FFA_Decode_faces = logPwrAll(:,:,1)';
            logPwr_PFC_Decode_faces = logPwrAll(:,:,2)';
            
            
            % DECODE GC
            GC_Decode_PFC_to_FFA_faces = GCAll(:,:,2)';
            GC_Decode_FFA_to_PFC_faces = GCAll(:,:,3)';
            
            logGC_Decode_PFC_to_FFA_faces = log_GCAll(:,:,1)'; %want in form of features (GC) x trials
            logGC_Decode_FFA_to_PFC_faces = log_GCAll(:,:,2)'; %want in form of features (GC) x trials
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %SAVE ALL THESE IN A STRUCTURE
            data.coh_faces = coh_Decode_faces;
            data.pwr_FFA_faces = pwr_FFA_Decode_faces;
            data.pwr_PFC_faces = pwr_PFC_Decode_faces;
            data.Logpwr_FFA_faces = logPwr_FFA_Decode_faces;
            data.Logpwr_PFC_faces = logPwr_PFC_Decode_faces;
            data.GC_PFC_to_FFA_faces = GC_Decode_PFC_to_FFA_faces;
            data.GC_FFA_to_PFC_faces = GC_Decode_FFA_to_PFC_faces;
            data.logGC_PFC_to_FFA_faces = logGC_Decode_PFC_to_FFA_faces;
            data.logGC_FFA_to_PFC_faces = logGC_Decode_FFA_to_PFC_faces;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Extract RANDOM coherence, power and GC for single-trial estimates %%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            load(['allrandom_CohPower_GC_SYNTHETIC_XX=' num2str(this_a) '_XY=' num2str(this_c) '_YY=' num2str(this_d) '.mat'])
            clear cohData cohAll pwrData pwrAll GCData GCAll
            %% COHERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %First look at coherence per trial
            freqs =  allGranger(1).freq
            freq_roi = 0 <= freqs & freqs < maxFreq
            freqWindow = freqs(freq_roi)
            freqIdxs = 1:length(freqWindow)
            
            
            for trialNo = 1:size(allCoherenceInfo,1)
                for thisFreq = 1:length(freqWindow)
                    cohData = allCoherenceInfo(trialNo).cohspctrm(:,:,thisFreq)
                    cohAll(trialNo,thisFreq,:) = cohData(:)'; % coherence (form trials x freq x cohVals)
                end
            end
            
            mean_Coh = squeeze(mean(cohAll,1))
            std_Coh = squeeze(std(cohAll))
            
            mean_Coh_b = mean_Coh(:,2)
            std_Coh_b = std_Coh(:,2)
            sem_Coh_b = std_Coh_b/(sqrt(length(std_Coh_b)))
            
            figure(1)
            errorbar(mean_Coh_b,sem_Coh_b, 'LineWidth', 4); hold on;
            legend({'FACES', 'RANDOM'});
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% POWER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for ss = 1:size(allPowerInfo,1)
                pwrData = allPowerInfo(ss).powspctrm(:,1:length(freqWindow))
                pwr_All(ss,:,:) = pwrData'; % power (form trials x freq x pwrPerCond)
            end
            
            mean_pwr = squeeze(mean(pwr_All,1))
            std_pwr = squeeze(std(pwr_All))
            sem_pwr = std_pwr/(sqrt(length(std_pwr)))
            
            %Log power
            logPwrAll = log10(pwr_All);
            mean_log_pwr = squeeze(mean(logPwrAll,1))
            std_log_pwr = squeeze(std(logPwrAll))
            sem_log_GC = std_log_pwr/(sqrt(length(std_log_pwr)))
            
            
            figure(21)
            errorbar(mean_pwr(:,1),sem_pwr(:,1), 'LineWidth', 4); hold on;
            legend({'FACES','RANDOM'});
            
            figure(22)
            errorbar(mean_pwr(:,2),sem_pwr(:,2), 'LineWidth', 4); hold on;
            legend({'FACES','RANDOM'});
            
            figure(26)
            imagesc(mean_pwr'); hold on;
            set(gca, 'FontSize', 30); xlabel('Frequency (Hz)','FontSize',30);
            ylabel('Power (mean across trials)','FontSize',30);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {freqWindow(1),freqWindow(10),freqWindow(20)}, 'YTick',[1,2],'YTickLabel',{'FFA','PFC'});
            xlim([1 length(mean_Coh_b)]); title('POWER - Left  FFA - RANDOM TRIALS', 'FontSize', 20); hold on;
            
            
            figure(24)
            errorbar(mean_log_pwr(:,1),sem_log_GC(:,1), 'LineWidth', 4); hold on;
            legend({'FACES','RANDOM'});
            
            figure(25)
            errorbar(mean_log_pwr(:,2),sem_log_GC(:,2), 'LineWidth', 4); hold on;
            legend({'FACES','RANDOM'});
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% GC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for trialNo = 1:size(allGranger,1)
                for thisFreq = 1:length(freqWindow)
                    GCData = allGranger(trialNo).grangerspctrm(:,:,thisFreq)
                    GCAll(trialNo,thisFreq,:) = GCData(:)'; % Granger (form trials x freq x GCVals)
                end
            end
            
            %Gather mean and std of GC
            mean_GC = squeeze(mean(GCAll,1))
            std_GC = squeeze(std(GCAll))
            sem_GC = std_GC/(sqrt(length(std_GC)))
            
            %Log GC
            useGC = GCAll(:,:,2:3)
            log_GCAll = log10(useGC)
            mean_log_GC = squeeze(mean(log_GCAll,1))
            std_log_GC = squeeze(std(log_GCAll))
            sem_log_GC = std_log_GC/(sqrt(length(std_log_GC)))
            
            %PLOT
            figure(31)
            errorbar(mean_GC(:,2),sem_GC(:,2), 'LineWidth', 4); hold on;
            legend({'FACES', 'RANDOM'});
            
            figure(32)
            errorbar(mean_GC(:,3),sem_GC(:,3), 'LineWidth', 4); hold on;
            legend({'FACES', 'RANDOM'});
            
            
            figure(33)
            errorbar(mean_log_GC(:,1),sem_log_GC(:,1), 'LineWidth', 4); hold on;
            legend({'FACES', 'RANDOM'});
            
            figure(34)
            errorbar(mean_log_GC(:,2),sem_log_GC(:,2), 'LineWidth', 4); hold on;
            legend({'FACES', 'RANDOM'});
            
            
            %Gather GC and logGC data for histogram plot
            plot_logGC_first_random = log_GCAll(:,5,2);
            plot_GC_first_random = GCAll(:,1,3);
            
            figure(8888)
            subplot(3,1,1)
            h1 = histogram(plot_logGC_first_faces(:), 'FaceColor', 'b'); hold on;
            title('Distribution of LOG for FFA to PFC at Freq#1'); set(gca, 'FontSize', 20);
            xlim([-5 5]); ylim([0 30]); h1.BinWidth = 0.2; axis off
            
            subplot(3,1,3)
            h2 = histogram(plot_logGC_first_random(:), 'FaceColor', 'g'); hold on;
            xlim([-5 5]); ylim([0 30]); h2.BinWidth = 0.2; hold off;
            legend({'FACES','RANDOM'}); axis off
            
            figure(88)
            h1 = histogram(plot_GC_first_faces(:), 'FaceColor', 'b'); hold on;
            title('Distribution of GC for FFA to PFC at Freq#1'); set(gca, 'FontSize', 20);
            xlim([0 10]); ylim([0 20]); h1.BinWidth = 0.2;
            h2 = histogram(plot_GC_first_random(:), 'FaceColor', 'g'); hold on;
            xlim([0 10]); ylim([0 20]); h2.BinWidth = 0.2; hold off;
            legend({'FACES','RANDOM'});
            
            %% Get data for decoding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % % DECODE COHERENCE
            coh_Decode_random = cohAll(:,:,3)'; % want in form of features (coh at freq) x trials
            
            % DECODE POWER
            pwr_FFA_Decode_random = pwr_All(:,:,1)'; % want in form of features(power) x trials
            pwr_PFC_Decode_random = pwr_All(:,:,2)'; % want in form of features(power) x trials
            
            logPwr_FFA_Decode_random = logPwrAll(:,:,1)';
            logPwr_PFC_Decode_random = logPwrAll(:,:,2)';
            
            % DECODE GC
            GC_Decode_PFC_to_FFA_random = GCAll(:,:,2)'; %want in form of features (GC) x trials
            GC_Decode_FFA_to_PFC_random = GCAll(:,:,3)'; %want in form of features (GC) x trials
            
            logGC_Decode_PFC_to_FFA_random = log_GCAll(:,:,1)'; %want in form of features (GC) x trials
            logGC_Decode_FFA_to_PFC_random = log_GCAll(:,:,2)'; %want in form of features (GC) x trials
            
            
            %SAVE ALL THESE IN A STRUCTURE
            data.coh_random = coh_Decode_random;
            data.pwr_FFA_random = pwr_FFA_Decode_random;
            data.pwr_PFC_random = pwr_PFC_Decode_random;
            data.Logpwr_FFA_random = logPwr_FFA_Decode_random;
            data.Logpwr_PFC_random = logPwr_PFC_Decode_random;
            data.GC_PFC_to_FFA_random = GC_Decode_PFC_to_FFA_random;
            data.GC_FFA_to_PFC_random = GC_Decode_FFA_to_PFC_random;
            data.logGC_PFC_to_FFA_random = logGC_Decode_PFC_to_FFA_random;
            data.logGC_FFA_to_PFC_random = logGC_Decode_FFA_to_PFC_random;
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% SINGLE FEATURE DECODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Decode COHERENCE usng libSVM
            nReps = 10;
            costRange = [-5] %-20:5:15;
            zscore = 1;
            
            [decodability_coh, cv_acc] = decode_singleFeature_libSVM(coh_Decode_faces, coh_Decode_random, nReps, costRange, zscore);
            
            %% Visualise
            decodability_mean = squeeze(mean(decodability_coh,2)); % mean decoding of each channel after 10 repetitions
            decodability_std = std(decodability_coh'); % std across each channel after repetitions
            
            results.singleCoh = decodability_coh;
            results.singleCoh_mean = decodability_mean;
            results.singleCoh_std = decodability_std;
            
            
            figure(11)
            errorbar(decodability_mean,decodability_std, 'LineWidth', 4); hold on;
            set(gca, 'FontSize', 30); xlabel('Frequency (Hz)','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',30);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {freqWindow(1),freqWindow(10),freqWindow(20)})
            xlim([1 length(decodability_mean)]);
            title('COHERENCE: Classification Accuracy - Face v Random Trials at different frequency bands', 'FontSize', 20);
            plot(1:length(decodability_mean)+1, ones(1,length(decodability_mean)+1)*50, '-r', 'LineWidth', 2); hold on;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Decode POWER usng libSVM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nReps = 10;
            costRange = [-5] %-20:5:15;
            zscore = 1;
            
            [decodability_coh, cv_acc] = decode_singleFeature_libSVM(pwr_FFA_Decode_faces, pwr_FFA_Decode_random, nReps, costRange, zscore);
            
            %% Visualise
            decodability_mean = squeeze(mean(decodability_coh,2)); % mean decoding of each channel after 10 repetitions
            decodability_std = std(decodability_coh'); % std across each channel after repetitions
            
            %SAVE THESE IN ONE STRUCTURE TOGETHER
            results.singlePwr_FFA = decodability_coh;
            results.singlePwr_FFA_mean = decodability_mean;
            results.singlePwr_FFA_std = decodability_std;
            
            
            figure(211)
            errorbar(decodability_mean,decodability_std, 'LineWidth', 4); hold on;
            set(gca, 'FontSize', 30); xlabel('Frequency (Hz)','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',30);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {freqWindow(1),freqWindow(10),freqWindow(20)})
            xlim([1 length(decodability_mean)]);
            title('POWER - L PFC: Classification Accuracy - Face v Random Trials at different frequency bands', 'FontSize', 20);
            plot(1:length(decodability_mean)+1, ones(1,length(decodability_mean)+1)*50, '-r', 'LineWidth', 2); hold on;
            
            [decodability_coh, cv_acc] = decode_singleFeature_libSVM(pwr_PFC_Decode_faces, pwr_PFC_Decode_random, nReps, costRange, zscore);
            
            %% Visualise
            decodability_mean = squeeze(mean(decodability_coh,2)); % mean decoding of each channel after 10 repetitions
            decodability_std = std(decodability_coh'); % std across each channel after repetitions
            
            figure(211)
            errorbar(decodability_mean,decodability_std, 'LineWidth', 4); hold on;
            legend({'FFA','PFC'});
            
            %SAVE THESE DATA
            results.singlePwr_PFC = decodability_coh;
            results.singlePwr_PFC_mean = decodability_mean;
            results.singlePwr_PFC_std = decodability_std;
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Decode log POWER usng libSVM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nReps = 10;
            costRange = [-5] %-20:5:15;
            zscore = 1;
            
            [decodability_coh, cv_acc] = decode_singleFeature_libSVM(logPwr_FFA_Decode_faces, logPwr_FFA_Decode_random, nReps, costRange, zscore);
            
            %% Visualise
            decodability_mean = squeeze(mean(decodability_coh,2)); % mean decoding of each channel after 10 repetitions
            decodability_std = std(decodability_coh'); % std across each channel after repetitions
            
            %SAVE THESE IN ONE STRUCTURE TOGETHER
            results.singleLogPwr_FFA = decodability_coh;
            results.singleLogPwr_FFA_mean = decodability_mean;
            results.singleLogPwr_FFA_std = decodability_std;
            
            figure(214)
            errorbar(decodability_mean,decodability_std, 'LineWidth', 4); hold on;
            set(gca, 'FontSize', 30); xlabel('Frequency (Hz)','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',30);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {freqWindow(1),freqWindow(10),freqWindow(20)})
            xlim([1 length(decodability_mean)]);
            title('LOG POWER - L PFC: Classification Accuracy - Face v Random Trials at different frequency bands', 'FontSize', 20);
            plot(1:length(decodability_mean)+1, ones(1,length(decodability_mean)+1)*50, '-r', 'LineWidth', 2); hold on;
            
            [decodability_coh, cv_acc] = decode_singleFeature_libSVM(logPwr_PFC_Decode_faces, logPwr_PFC_Decode_random, nReps, costRange, zscore);
            
            %% Visualise
            decodability_mean = squeeze(mean(decodability_coh,2)); % mean decoding of each channel after 10 repetitions
            decodability_std = std(decodability_coh'); % std across each channel after repetitions
            
            figure(214)
            errorbar(decodability_mean,decodability_std, 'LineWidth', 4); hold on;
            legend({'FFA','PFC'});
            
            %SAVE THESE DATA
            results.singleLogPwr_PFC = decodability_coh;
            results.singleLogPwr_PFC_mean = decodability_mean;
            results.singleLogPwr_PFC_std = decodability_std;
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Decode GC usng libSVM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nReps = 10;
            costRange = [-5] %-20:5:15;
            zscore = 1;
            
            [decodability_coh, cv_acc] = decode_singleFeature_libSVM(GC_Decode_FFA_to_PFC_faces, GC_Decode_FFA_to_PFC_random, nReps, costRange, zscore);
            
            %% Visualise
            decodability_mean = squeeze(mean(decodability_coh,2)); % mean decoding of each channel after 10 repetitions
            decodability_std = std(decodability_coh'); % std across each channel after repetitions
            
            %SAVE THESE IN ONE STRUCTURE TOGETHER
            results.singleGC_FFA_to_PFC = decodability_coh;
            results.singleGC_FFA_to_PFC_mean = decodability_mean;
            results.singleGC_FFA_to_PFC_std = decodability_std;
            
            
            figure(311)
            errorbar(decodability_mean,decodability_std, 'LineWidth', 4); hold on;
            set(gca, 'FontSize', 30); xlabel('Frequency (Hz)','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',30);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {freqWindow(1),freqWindow(10),freqWindow(20)})
            xlim([1 length(decodability_mean)]);
            title('GC - L FFA and L PFC: Classification Accuracy - Face v Random Trials at different frequency bands', 'FontSize', 20);
            plot(1:length(decodability_mean)+1, ones(1,length(decodability_mean)+1)*50, '-r', 'LineWidth', 2); hold on;
            
            
            [decodability_coh, cv_acc] = decode_singleFeature_libSVM(GC_Decode_PFC_to_FFA_faces, GC_Decode_PFC_to_FFA_random, nReps, costRange, zscore);
            
            decodability_mean = squeeze(mean(decodability_coh,2)); % mean decoding of each channel after 10 repetitions
            decodability_std = std(decodability_coh'); % std across each channel after repetitions
            
            figure(311)
            errorbar(decodability_mean,decodability_std, 'LineWidth', 4); hold on;
            legend({'FFA to PFC','','PFC to FFA'});
            
            %SAVE THESE DATA
            
            results.singleGC_PFC_to_FFA = decodability_coh;
            results.singleGC_PFC_to_FFA_mean = decodability_mean;
            results.singleGC_PFC_to_FFA_std = decodability_std;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Decode log GC usng libSVM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nReps = 10;
            costRange = [-5] %-20:5:15;
            zscore = 1;
            
            [decodability_GC_FFA_PFC, cv_acc] = decode_singleFeature_libSVM(logGC_Decode_FFA_to_PFC_faces, logGC_Decode_FFA_to_PFC_random, nReps, costRange, zscore);
            
            %% Visualise
            decodability_mean = squeeze(mean(decodability_GC_FFA_PFC,2)); % mean decoding of each channel after 10 repetitions
            decodability_std = std(decodability_GC_FFA_PFC'); % std across each channel after repetitions
            
            
            %SAVE THESE IN ONE STRUCTURE TOGETHER
            results.singleLogGC_FFA_to_PFC = decodability_GC_FFA_PFC;
            results.singleLogGC_FFA_to_PFC_mean = decodability_mean;
            results.singleLogGC_FFA_to_PFC_std = decodability_std;
            
            figure(3111)
            errorbar(decodability_mean,decodability_std, 'LineWidth', 4); hold on;
            set(gca, 'FontSize', 30); xlabel('Frequency (Hz)','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',30);
            set(gca, 'XTick', [1,10,20], 'XTickLabel', {freqWindow(1),freqWindow(10),freqWindow(20)})
            xlim([1 length(decodability_mean)]);
            title('Log GC - L FFA and L PFC: Classification Accuracy - Face v Random Trials at different frequency bands', 'FontSize', 20);
            plot(1:length(decodability_mean)+1, ones(1,length(decodability_mean)+1)*50, '-r', 'LineWidth', 2); hold on;
            
            
            [decodability_GC_PFC_FFA, cv_acc] = decode_singleFeature_libSVM(logGC_Decode_PFC_to_FFA_faces, logGC_Decode_PFC_to_FFA_random, nReps, costRange, zscore);
            
            decodability_mean = squeeze(mean(decodability_GC_PFC_FFA,2)); % mean decoding of each channel after 10 repetitions
            decodability_std = std(decodability_GC_PFC_FFA'); % std across each channel after repetitions
            
            results.singleLogGC_PFC_to_FFA = decodability_GC_PFC_FFA;
            results.singleLogGC_PFC_to_FFA_mean = decodability_mean;
            results.singleLogGC_PFC_to_FFA_std = decodability_std;
            
            
            figure(3111)
            errorbar(decodability_mean,decodability_std, 'LineWidth', 4); hold on;
            legend({'FFA to PFC','','PFC to FFA'});
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% MULTI FEATURE DECODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Decode COHERENCE usng libSVM
            nReps = 10;
            costRange = [-5] %-20:5:15;
            zscore = 1;
            
            [decodability_coh, cv_acc] = decode_multiFeature_libSVM(coh_Decode_faces, coh_Decode_random, nReps, costRange, zscore);
            
            %% Visualise
            decodability_mean = squeeze(mean(decodability_coh)); % mean decoding of each channel after 10 repetitions
            decodability_std = std(decodability_coh'); % std across each channel after repetitions
            
            results.multiCoh = decodability_coh;
            results.multiCoh_mean = decodability_mean;
            results.multiCoh_std = decodability_std;
            
            figure(1111)
            bar(decodability_mean); hold on;
            ylim([45 100]); xlim([0 length(decodability_mean)+1])
            errorbar(1:length(decodability_mean),decodability_mean, decodability_std,decodability_std)
            set(gca, 'FontSize', 30); xlabel('Coherence (multifeature)','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',20);
            plot(0:length(decodability_mean)+1, ones(1,length(decodability_mean)+2)*50, '-r', 'LineWidth', 2)
            title('COHERENCE Classification Accuracy - Face v Random Trials - Multifeature', 'FontSize', 20);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Decode POWER usng libSVM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nReps = 10;
            costRange = [-5] %-20:5:15;
            zscore = 1;
            
            [decodability_logPwr_FFA, cv_acc] = decode_multiFeature_libSVM(pwr_FFA_Decode_faces, pwr_FFA_Decode_random, nReps, costRange, zscore);
            
            %% Visualise
            decodability_mean_FFA = squeeze(mean(decodability_logPwr_FFA)); % mean decoding of each channel after 10 repetitions
            decodability_std_FFA = std(decodability_logPwr_FFA'); % std across each channel after repetitions
            
            [decodability_logPwr_PFC, cv_acc] = decode_multiFeature_libSVM(pwr_PFC_Decode_faces, pwr_PFC_Decode_random, nReps, costRange, zscore);
            
            %% Visualise
            decodability_mean_PFC = squeeze(mean(decodability_logPwr_PFC)); % mean decoding of each channel after 10 repetitions
            decodability_std_PFC = std(decodability_logPwr_PFC'); % std across each channel after repetitions
            
            decode_logpwr_FFA_and_PFC = [decodability_mean_FFA, decodability_mean_PFC]
            decode_logpwr_FFA_PFC_std = [decodability_std_FFA, decodability_std_PFC]
            
            %SAVE THESE IN ONE STRUCTURE TOGETHER
            results.multiPwr_FFA = decodability_logPwr_FFA;
            results.multiPwr_FFA_mean = decodability_mean_FFA;
            results.multiPwr_FFA_std = decodability_std_FFA;
            
            results.multiPwr_PFC = decodability_logPwr_PFC;
            results.multiPwr_PFC_mean = decodability_mean_PFC;
            results.multiPwr_PFC_std = decodability_std_PFC;
            
            
            figure(2222)
            bar(decode_logpwr_FFA_and_PFC); hold on;
            ylim([45 70]); xlim([0 length(decode_logpwr_FFA_and_PFC)+1])
            errorbar(1:length(decode_logpwr_FFA_and_PFC),decode_logpwr_FFA_and_PFC, decode_logpwr_FFA_PFC_std,decode_logpwr_FFA_PFC_std, 'LineWidth', 4)
            set(gca, 'FontSize', 30); xlabel('Power (multifeature)','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',20);
            plot(0:length(decode_logpwr_FFA_PFC_std)+1, ones(1,length(decode_logpwr_FFA_PFC_std)+2)*50, '-r', 'LineWidth', 2); hold on;
            set(gca, 'XTick', [1, 2], 'XTickLabel',{'FFA','PFC'});
            title('POWER Classification Accuracy - Face v Random Trials - Multifeature', 'FontSize', 20);
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Decode log POWER usng libSVM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nReps = 10;
            costRange = [-5] %-20:5:15;
            zscore = 1;
            
            [decodability_logPwr_FFA, cv_acc] = decode_multiFeature_libSVM(logPwr_FFA_Decode_faces, logPwr_FFA_Decode_random, nReps, costRange, zscore);
            
            %% Visualise
            decodability_mean_FFA = squeeze(mean(decodability_logPwr_FFA)); % mean decoding of each channel after 10 repetitions
            decodability_std_FFA = std(decodability_logPwr_FFA'); % std across each channel after repetitions
            
            [decodability_logPwr_PFC, cv_acc] = decode_multiFeature_libSVM(logPwr_PFC_Decode_faces, logPwr_PFC_Decode_random, nReps, costRange, zscore);
            
            %% Visualise
            decodability_mean_PFC = squeeze(mean(decodability_logPwr_PFC)); % mean decoding of each channel after 10 repetitions
            decodability_std_PFC = std(decodability_logPwr_PFC'); % std across each channel after repetitions
            
            decode_logpwr_FFA_and_PFC = [decodability_mean_FFA, decodability_mean_PFC]
            decode_logpwr_FFA_PFC_std = [decodability_std_FFA, decodability_std_PFC]
            
            %SAVE THESE IN ONE STRUCTURE TOGETHER
            results.multiLogPwr_FFA = decodability_logPwr_FFA;
            results.multiLogPwr_FFA_mean = decodability_mean_FFA;
            results.multiLogPwr_FFA_std = decodability_std_FFA;
            
            results.multiLogPwr_PFC = decodability_logPwr_PFC;
            results.multiLogPwr_PFC_mean = decodability_mean_PFC;
            results.multiLogPwr_PFC_std = decodability_std_PFC;
            
            
            figure(21228)
            bar(decode_logpwr_FFA_and_PFC); hold on;
            ylim([45 100]); xlim([0 length(decode_logpwr_FFA_and_PFC)+1])
            errorbar(1:length(decode_logpwr_FFA_and_PFC),decode_logpwr_FFA_and_PFC, decode_logpwr_FFA_PFC_std,decode_logpwr_FFA_PFC_std, 'LineWidth', 6)
            set(gca, 'FontSize', 30); xlabel('Power (multifeature)','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',20);
            plot(0:length(decode_logpwr_FFA_PFC_std)+1, ones(1,length(decode_logpwr_FFA_PFC_std)+2)*50, '-r', 'LineWidth', 4); hold on;
            set(gca, 'XTick', [1, 2], 'XTickLabel',{'FFA','PFC'});
            %    title('LOG POWER Classification Accuracy - Face v Random Trials - Multifeature', 'FontSize', 20);
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Decode GC usng libSVM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            nReps = 10;
            costRange = [-5] %-20:5:15;
            zscore = 1;
            
            [decodability_GC_FFA_to_PFC, cv_acc] = decode_multiFeature_libSVM(GC_Decode_FFA_to_PFC_faces, GC_Decode_FFA_to_PFC_random, nReps, costRange, zscore);
            
            %% Visualise
            decodability_mean_FFA_PFC = squeeze(mean(decodability_GC_FFA_to_PFC)); % mean decoding of each channel after 10 repetitions
            decodability_std_FFA_PFC = std(decodability_GC_FFA_to_PFC'); % std across each channel after repetitions
            
            
            [decodability_GC_PFC_to_FFA, cv_acc] = decode_multiFeature_libSVM(GC_Decode_PFC_to_FFA_faces, GC_Decode_PFC_to_FFA_random, nReps, costRange, zscore);
            
            decodability_mean_PFC_FFA = squeeze(mean(decodability_GC_PFC_to_FFA)); % mean decoding of each channel after 10 repetitions
            decodability_std_PFC_FFA = std(decodability_GC_PFC_to_FFA'); % std across each channel after repetitions
            
            
            decode_GC_FFA_PFC = [decodability_mean_FFA_PFC, decodability_mean_PFC_FFA]
            decode_GC_FFA_PFC_std = [decodability_std_FFA_PFC, decodability_std_PFC_FFA]
            
            
            %SAVE THESE IN ONE STRUCTURE TOGETHER
            results.multiGC_FFA_to_PFC = decodability_GC_FFA_to_PFC;
            results.multiGC_FFA_to_PFC_mean = decodability_mean_FFA_PFC;
            results.multiGC_FFA_to_PFC_std = decodability_std_FFA_PFC;
            
            results.multiGC_PFC_to_FFA = decodability_GC_PFC_to_FFA;
            results.multiGC_PFC_to_FFA_mean = decodability_mean_PFC_FFA;
            results.multiGC_PFC_to_FFA_std = decodability_std_PFC_FFA;
            
            
            figure(3333)
            bar(decode_GC_FFA_PFC); hold on;
            ylim([45 100]); xlim([0 length(decode_GC_FFA_PFC)+1])
            errorbar(1:length(decode_GC_FFA_PFC),decode_GC_FFA_PFC, decode_GC_FFA_PFC_std,decode_GC_FFA_PFC_std)
            set(gca, 'FontSize', 30); xlabel('GC (multifeature)','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',20);
            plot(0:length(decode_GC_FFA_PFC_std)+1, ones(1,length(decode_GC_FFA_PFC_std)+2)*50, '-r', 'LineWidth', 2); hold on;
            set(gca, 'XTick', [1, 2], 'XTickLabel',{'FFA to PFC','PFC to FFA'});
            title('GC Classification Accuracy - Face v Random Trials - Multifeature', 'FontSize', 20);
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             %% Decode log GC usng libSVM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %             nReps = 10;
            %             costRange = [-5] %-20:5:15;
            %             zscore = 1;
            %
            %             [decodability_FFA_PFC, cv_acc] = decode_multiFeature_libSVM(logGC_Decode_FFA_to_PFC_faces, logGC_Decode_FFA_to_PFC_random, nReps, costRange, zscore);
            %
            %             %% Visualise
            %             decodability_mean_FFA_PFC = squeeze(mean(decodability_FFA_PFC)); % mean decoding of each channel after 10 repetitions
            %             decodability_std_FFA_PFC = std(decodability_FFA_PFC'); % std across each channel after repetitions
            %
            %
            %             [decodability_PFC_FFA, cv_acc] = decode_multiFeature_libSVM(logGC_Decode_PFC_to_FFA_faces, logGC_Decode_PFC_to_FFA_random, nReps, costRange, zscore);
            %
            %             decodability_mean_PFC_FFA = squeeze(mean(decodability_PFC_FFA)); % mean decoding of each channel after 10 repetitions
            %             decodability_std_PFC_FFA = std(decodability_PFC_FFA'); % std across each channel after repetitions
            %
            %
            %             decode_GC_FFA_PFC = [decodability_mean_FFA_PFC, decodability_mean_PFC_FFA]
            %             decode_GC_FFA_PFC_std = [decodability_std_FFA_PFC, decodability_std_PFC_FFA]
            %
            %
            %             %SAVE THESE IN ONE STRUCTURE TOGETHER
            %             results.multiLogGC_FFA_to_PFC = decodability_FFA_PFC;
            %             results.multiLogGC_FFA_to_PFC_mean = decodability_mean_FFA_PFC;
            %             results.multiLogGC_FFA_to_PFC_std = decodability_std_FFA_PFC;
            %
            %             results.multiLogGC_PFC_to_FFA = decodability_PFC_FFA;
            %             results.multiLogGC_PFC_to_FFA_mean = decodability_mean_PFC_FFA;
            %             results.multiLogGC_PFC_to_FFA_std = decodability_std_PFC_FFA;
            %
            %
            %             %PLOT
            %             figure(33333)
            %             bar(decode_GC_FFA_PFC); hold on;
            %             ylim([45 100]); xlim([0 length(decode_GC_FFA_PFC)+1])
            %             errorbar(1:length(decode_GC_FFA_PFC),decode_GC_FFA_PFC, decode_GC_FFA_PFC_std,decode_GC_FFA_PFC_std)
            %             set(gca, 'FontSize', 30); xlabel('Log GC (multifeature)','FontSize',30); ylabel('Classification Accuracy (%)','FontSize',20);
            %             plot(0:length(decode_GC_FFA_PFC_std)+1, ones(1,length(decode_GC_FFA_PFC_std)+2)*50, '-r', 'LineWidth', 2); hold on;
            %             set(gca, 'XTick', [1, 2], 'XTickLabel',{'FFA to PFC','PFC to FFA'});
            %             title('Log GC Classification Accuracy - Face v Random Trials - Multifeature', 'FontSize', 20);
            
            %% SAVE DATA AND RESULTS
            save(savename, 'results','data')
            %   clear results; clear data;
            %   close all;
        end
    end
end

%PLOT ALL IN ONE FIGURE
%Source waveforms
mD_sWave_X = mean(decodability_sWaves_nodeX)
stdD_sWave_X = std(decodability_sWaves_nodeX)
mD_sWave_Y = mean(decodability_sWaves_nodeY)
stdD_sWave_Y = std(decodability_sWaves_nodeY)

all_sWave_X_Y = [mD_sWave_X,mD_sWave_Y]
all_sWave_std_X_Y = [stdD_sWave_X, stdD_sWave_Y]+0.0001;

% log power
mD_logPwr_PFC = mean(decodability_logPwr_PFC)
stdD_logPwr_PFC = std(decodability_logPwr_PFC)
mD_logPwr_FFA = mean(decodability_logPwr_FFA)
stdD_logPwr_FFA = std(decodability_logPwr_FFA)

all_logPwr_FFA_PFC = [mD_logPwr_PFC, mD_logPwr_FFA]
all_logPwr_std_FFA_PFC = [stdD_logPwr_PFC, stdD_logPwr_FFA]

%coherence
mD_coh = mean(decodability_coh)
stdD_coh = std(decodability_coh)

%GC to PFC
mD_GC_FG_PFC = mean(decodability_GC_FFA_to_PFC)
stdD_GC_FG_PFC = std(decodability_GC_FFA_to_PFC)

%GC from PFC
mD_GC_PFC_FG = mean(decodability_GC_PFC_to_FFA)
stdD_GC_PFC_FG = std(decodability_GC_PFC_to_FFA)

% PLOT (bar plot)
figure(343)
subplot(2,3,1)
bar(all_sWave_X_Y); hold on;
ylim([45 100]); xlim([0 length(all_sWave_X_Y)+1])
errorbar(1:length(all_sWave_X_Y),all_sWave_X_Y, all_sWave_std_X_Y)
set(gca, 'FontSize', 30); ylabel('Classification Accuracy (%)','FontSize',20);
plot(0:length(decode_GC_FFA_PFC_std)+1, ones(1,length(decode_GC_FFA_PFC_std)+2)*50, '-r', 'LineWidth', 2); hold on;
set(gca, 'XTick', [1, 2], 'XTickLabel',{'X','Y'});

subplot(2,3,2)
bar(all_logPwr_FFA_PFC); hold on;
ylim([45 100]); xlim([0 length(all_logPwr_FFA_PFC)+1])
errorbar(1:length(all_logPwr_FFA_PFC),all_logPwr_FFA_PFC, all_logPwr_std_FFA_PFC, 'LineWidth',3)
%set(gca, 'FontSize', 30); ylabel('Classification Accuracy (%)','FontSize',20);
plot(0:length(decode_GC_FFA_PFC_std)+1, ones(1,length(decode_GC_FFA_PFC_std)+2)*50, '-r', 'LineWidth', 2); hold on;
set(gca, 'XTick', [1, 2], 'XTickLabel',{'X','Y'},'FontSize',30 ,'YTickLabel',[]);

subplot(2,3,3)
bar(mD_coh); hold on;
ylim([45 100]); xlim([0 length(mD_coh)+1])
errorbar(1:length(mD_coh),mD_coh, stdD_coh, 'LineWidth',3)
%set(gca, 'FontSize', 30); ylabel('Classification Accuracy (%)','FontSize',20);
plot(0:length(decode_GC_FFA_PFC_std)+1, ones(1,length(decode_GC_FFA_PFC_std)+2)*50, '-r', 'LineWidth', 2); hold on;
set(gca, 'XTick', [1], 'XTickLabel',{'X&Y'},'FontSize',30 ,'YTickLabel',[]);

subplot(2,3,4)
bar(mD_GC_FG_PFC); hold on;
ylim([45 100]); xlim([0 length(mD_GC_FG_PFC)+1])
errorbar(1:length(mD_GC_FG_PFC),mD_GC_FG_PFC, stdD_GC_FG_PFC, 'LineWidth',3)
%set(gca, 'FontSize', 30); ylabel('Classification Accuracy (%)','FontSize',20);
plot(0:length(decode_GC_FFA_PFC_std)+1, ones(1,length(decode_GC_FFA_PFC_std)+2)*50, '-r', 'LineWidth', 2); hold on;
set(gca, 'XTick', [1], 'XTickLabel',{'X-->Y'},'FontSize',30 ,'YTickLabel',[]);

subplot(2,3,5)
bar(mD_GC_PFC_FG); hold on;
ylim([45 100]); xlim([0 length(mD_GC_PFC_FG)+1])
errorbar(1:length(mD_GC_PFC_FG),mD_GC_PFC_FG, stdD_GC_PFC_FG, 'LineWidth',3)
%set(gca, 'FontSize', 30); ylabel('Classification Accuracy (%)');
plot(0:length(decode_GC_FFA_PFC_std)+1, ones(1,length(decode_GC_FFA_PFC_std)+2)*50, '-r', 'LineWidth', 2); hold on;
set(gca, 'XTick', [1], 'XTickLabel',{'Y-->X'},'FontSize',30 ,'YTickLabel',[]);
