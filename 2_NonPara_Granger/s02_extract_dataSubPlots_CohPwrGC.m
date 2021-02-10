%% Visualise and decode NP coherence, power and Granger Causality results

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

close all; clear all;
counterL = 1;
counterP = 1;
addpath(genpath('/Users/egrow1/Desktop/NEW_Redo_OriginalChans_SHAFTO/decode_libSVM'))
%Add the libSVM decoding library to the path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT: Settings and filenames for extracting ROI data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filepath = ['/Users/egrow1/Desktop/ALL_SW_noF_0_to_500ms_EBR/']; %set file path
%filepathList = ['/Users/egrow1/Desktop/NEW_Redo_OriginalChans_SHAFTO/SourceResults_MASSIVE/'];
filePathCoords = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/1_EEG_PreProcessing/'];

% SELECT THE TWO NODES TO DECODE
possLowLevlCoords = {'L_FFA'}%, 'R_FFA','L_Occ','R_Occ'} %'main' PFC coords
possPFCCoords = {'L_PFC'}%,'R_PFC','L_Inf_PFC','R_Inf_PFC'}; %,'LeftPar','LeftOcc', ...
maxFreq = 55; % maximum frequency to include (i.e. 0 to maxFreq Hz)
countA = 1; countB = 1; countBB = 1; countC = 1; %setup indexes

for lowIdx = 1:length(possLowLevlCoords)
    
    useLowLvlCoord = possLowLevlCoords{lowIdx}
    
    for PFCIdx = 1:length(possPFCCoords)
        
        usePFCCoord = possPFCCoords{PFCIdx}
        
        % ASSIGN SAVE AND LOAD NAMES
        loadNameFaces = ['all_NoF_faces_CohPower_GC_' num2str(usePFCCoord) ...
            '_and_' num2str(useLowLvlCoord) '_0to500ms_5Tapsmo_TopVoxel_FINAL.mat']; %savename
        
        loadNameRandom = ['all_NoF_random_CohPower_GC_' num2str(usePFCCoord) ...
            '_and_' num2str(useLowLvlCoord) '_0to500ms_5Tapsmo_TopVoxel_FINAL.mat']; %savename
        
        savename = ['results_noF_Decode_CohPowerGC_' num2str(usePFCCoord) ...
            '_and_' num2str(useLowLvlCoord) '_251Hz_0to500ms_5Tapsmo_TopVoxel_FINAL.mat']; %savename
        
        % ASSIGN DECODING PARAMETERS (libSVM)
        nReps = 10; % number of repetitions of decoding steps
        costRange = [-5]; %cost range = 1 value (%-20:5:15) set to range if doing optimisation
        zscore = 1; % zscore decoding data (1 = 'yes', 0 = 'no')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ANALYSIS 1a: COHERENCE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cohAll_faces = [];  cohAll_random = [];
        %% Load Face Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load(loadNameFaces); %load face data
        
        %First look at coherence per trial
        freq_roi =  0 <= allGranger(1).freq & allGranger(1).freq < maxFreq
        freqWindow = allGranger(1).freq(freq_roi)
        
        for trialNo = 1:size(allCoherenceInfo,1)
            for thisFreq = 1:length(freqWindow)
                cohData_faces = allCoherenceInfo(trialNo).cohspctrm(:,:,thisFreq); %extract data
                cohAll_faces(trialNo,thisFreq,:) = cohData_faces(:)'; % coherence (form trials x freq x cohVals)
            end
        end
        
        cohAll_faces = cohAll_faces(:,:,2); % extract only value that matters
        
        %Find the mean coherence value across all trials (per frequency)
        mCoh_faces = squeeze(mean(cohAll_faces,1))'; %mean all coh vals
        stdCoh_faces = squeeze(std(cohAll_faces,1))'; %std all coh vals
        semCoh_faces = stdCoh_faces/(sqrt(length(cohAll_faces))); % SEM coherence between ROIs
        
        %PLOT COHERENCE VALUES: Faces
        figure(1)
        %thisPlot = usePlotIdxs(countA);
        subplot(4,4,countA)
        errorbar(mCoh_faces,semCoh_faces, 'LineWidth', 3); hold on; ylim([0.35 0.65]); xlim([1 length(mCoh_faces)]);
        if PFCIdx == 1 && lowIdx == 1
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),length(freqWindow)], 'XTickLabel', ...
                {ceil(freqWindow(1)),ceil(freqWindow(ceil(length(freqWindow)/2))),ceil(freqWindow(end))});
        else
            set(gca,'XTick',[]); set(gca,'YTick',[])
        end %set(gca, 'FontSize', 1); end
        %    title(['COH: ' num2str(useLowLvlCoord) ' and ' num2str(usePFCCoord)], 'FontSize', 15); %legend({'', 'TEST', '', 'TRAIN'});
        hold on;
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Load Random Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load(loadNameRandom); %load random data
        
        for trialNo = 1:size(allCoherenceInfo,1)
            for thisFreq = 1:length(freqWindow)
                cohData_random = allCoherenceInfo(trialNo).cohspctrm(:,:,thisFreq)
                cohAll_random(trialNo,thisFreq,:) = cohData_random(:)'; % coherence (form trials x freq x cohVals)
            end
        end
        
        cohAll_random = cohAll_random(:,:,2); % extract only value that matters
        
        %Find the mean coherence value across all trials (per frequency)
        mCoh_random = squeeze(mean(cohAll_random,1))'; %mean all coh vals
        stdCoh_random = squeeze(std(cohAll_random,1))'; %std all coh vals
        semCoh_random = stdCoh_random/(sqrt(length(cohAll_random))); % SEM coherence between ROIs
        
        figure(1)
        subplot(4,4,countA)
        errorbar(mCoh_random,semCoh_random, 'LineWidth', 3); hold on;
        
        
        %Coh data for face and random
        all_m_Coh = [mCoh_faces'; mCoh_random']
        all_sem_Coh = [semCoh_faces'; semCoh_random']
        
        figure(11111)
        subplot(4,4,countA)
        mseb((repmat(1:length(all_m_Coh),2,1)),all_m_Coh,all_sem_Coh,[],1); box on;
        
        if PFCIdx == 1 && lowIdx == 1
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15); ylabel('Coherence', 'FontSize', 15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),length(freqWindow)], 'XTickLabel', ...
                {ceil(freqWindow(1)),ceil(freqWindow(ceil(length(freqWindow)/2))),ceil(freqWindow(end))});
            xlim([1 length(all_m_Coh)]); ylim([0.35 0.65]);
        else
            set(gca,'XTick',[]); set(gca,'YTick',[]);  xlim([1 length(all_m_Coh)]);  ylim([0.35 0.65]);
        end
        countA = countA + 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ANALYSIS 3a: GRANGER CAUSALITY RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Load Face Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load(loadNameFaces); %load face data
        
        for trialNo = 1:size(allGranger,1)
            for thisFreq = 1:length(freqWindow)
                GCData_faces = allGranger(trialNo).grangerspctrm(:,:,thisFreq); % GC data
                GCAll_faces(trialNo,thisFreq,:) = GCData_faces(:)'; % Granger (form trials x freq x GCVals)
            end
        end
        
        %Gather mean and std of GC
        mGC_faces = squeeze(mean(GCAll_faces,1)); %MEAN GC over trials (at each frequency band)
        stdGC_faces = squeeze(std(GCAll_faces)); %STD GC over trials
        semGC_faces = stdGC_faces/(sqrt(length(GCAll_faces))); %SEM GC over trials
        
        %Log GC
        useGC_faces = GCAll_faces(:,:,2:3); %only examine directionaly GC coefficients (i.e. column 2 and 3)
        logGCAll_faces = log10(useGC_faces); % log transform the GC values (improve distribution)
        mLogGC_faces = squeeze(mean(logGCAll_faces,1)); %MEAN log GC values
        stdLogGC_faces = squeeze(std(logGCAll_faces)); % STD log GC values
        semLogGC_faces = stdLogGC_faces/(sqrt(length(logGCAll_faces))); % SEM log GC values
        
        GC_TO_PFC_faces = mGC_faces(:,3)
        GC_sem_TO_PFC_faces = semGC_faces(:,3)
        
        GC_FROM_PFC_faces = mGC_faces(:,2)
        GC_sem_FROM_PFC_faces = semGC_faces(:,2)
        
        figure(333)
        %thisPlot = usePlotIdxs(countA);
        subplot(4,4,countC)
        errorbar(GC_TO_PFC_faces,GC_sem_TO_PFC_faces, 'LineWidth', 3); hold on; ylim([0.1 0.5]); xlim([1 length(mCoh_faces)]);
        if PFCIdx == 1 && lowIdx == 1
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),length(freqWindow)], 'XTickLabel', ...
                {ceil(freqWindow(1)),ceil(freqWindow(ceil(length(freqWindow)/2))),ceil(freqWindow(end))});
        else set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]); xlim([1 length(mCoh_faces)]);
        end
        %  title(['GC to ' num2str(usePFCCoord) ' from:' num2str(useLowLvlCoord)], 'FontSize', 15); %legend({'', 'TEST', '', 'TRAIN'});
        hold on;
        
        figure(444)
        subplot(4,4,countC)
        errorbar(GC_FROM_PFC_faces,GC_sem_FROM_PFC_faces, 'LineWidth', 3); hold on; ylim([0.1 0.5]); xlim([1 length(mCoh_faces)]);
        if PFCIdx == 1 && lowIdx == 1
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),length(freqWindow)], 'XTickLabel', ...
                {ceil(freqWindow(1)),ceil(freqWindow(ceil(length(freqWindow)/2))),ceil(freqWindow(end))});
        else set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]); xlim([1 length(mCoh_faces)]);
        end
        %  title(['GC from ' num2str(usePFCCoord) ' to: ' num2str(useLowLvlCoord)], 'FontSize', 15); %legend({'', 'TEST', '', 'TRAIN'});
        hold on;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Load Random Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load(loadNameRandom); %load random data
        
        for trialNo = 1:size(allGranger,1)
            for thisFreq = 1:length(freqWindow)
                GCData_random = allGranger(trialNo).grangerspctrm(:,:,thisFreq); % extract GC data
                GCAll_random(trialNo,thisFreq,:) = GCData_random(:)'; % Granger (form trials x freq x GCVals)
            end
        end
        
        %Gather mean and std of GC
        mGC_random = squeeze(mean(GCAll_random,1))
        stdGC_random = squeeze(std(GCAll_random))
        semGC_random = stdGC_random/(sqrt(length(GCAll_random)))
        
        %Calculate Log GC
        useGC_random = GCAll_random(:,:,2:3)
        logGCAll_random = log10(useGC_random)
        mLogGC_random = squeeze(mean(logGCAll_random,1))
        stdLogGC_random = squeeze(std(logGCAll_random))
        semLogGC_random = stdLogGC_random/(sqrt(length(GCAll_random)))
        
        
        GC_TO_PFC_random = mGC_random(:,3)
        GC_sem_TO_PFC_random = semGC_random(:,3)
        
        GC_FROM_PFC_random = mGC_random(:,2)
        GC_sem_FROM_PFC_random = semGC_random(:,2)
        
        figure(333)
        subplot(4,4,countC)
        errorbar(GC_TO_PFC_random,GC_sem_TO_PFC_random, 'LineWidth', 3); hold on;
        
        figure(444)
        subplot(4,4,countC)
        errorbar(GC_FROM_PFC_random,GC_sem_FROM_PFC_random, 'LineWidth', 3); hold on;
        
        %Gather both data
        all_GC_to_PFC = [GC_TO_PFC_faces'; GC_TO_PFC_random'];
        all_GC_to_PFC_sem = [GC_sem_TO_PFC_faces'; GC_sem_TO_PFC_random'];
        
        all_GC_from_PFC = [GC_FROM_PFC_faces'; GC_FROM_PFC_random'];
        all_GC_from_PFC_sem = [GC_sem_FROM_PFC_faces'; GC_sem_FROM_PFC_random'];
        
        figure(31111)
        subplot(4,4,countC)
        mseb((repmat(1:length(all_GC_to_PFC),2,1)),all_GC_to_PFC,all_GC_to_PFC_sem,[],1); box on;
        if PFCIdx == 1 && lowIdx == 1
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),length(freqWindow)], 'XTickLabel', ...
                {ceil(freqWindow(1)),ceil(freqWindow(ceil(length(freqWindow)/2))),ceil(freqWindow(end))});
            xlim([1 length(all_m_Coh)]);  ylim([0.1 0.45]);
        else set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]); xlim([1 length(all_m_Coh)]); ylim([0.1 0.45]);
        end
        
        figure(41111)
        subplot(4,4,countC)
        mseb((repmat(1:length(all_GC_from_PFC),2,1)),all_GC_from_PFC,all_GC_from_PFC_sem,[],1);box on;
        if PFCIdx == 1 && lowIdx == 1
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),length(freqWindow)], 'XTickLabel', ...
                {ceil(freqWindow(1)),ceil(freqWindow(ceil(length(freqWindow)/2))),ceil(freqWindow(end))});
            xlim([1 length(all_m_Coh)]); ylim([0.1 0.45]);
        else set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]); xlim([1 length(all_m_Coh)]);  ylim([0.1 0.45]);
        end
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
        
        mPwr_faces = squeeze(mean(pwrAll_faces,1)); %MEAN power across trials
        stdPwr_faces = squeeze(std(pwrAll_faces)); % STDEV power across trials
        semPwr_faces = stdPwr_faces/(sqrt(length(pwrAll_faces))); %SEM power across trials
        
        
        %Log power
        logPwrAll_faces = log10(pwrAll_faces*10); %take log10 of power (improve distribution)
        mLogPwr_faces = squeeze(mean(logPwrAll_faces,1)); %MEAN of log transformed power
        stdLogPwr_faces = squeeze(std(logPwrAll_faces)); %STDEV of log transformed power
        semLogPwr_faces = stdLogPwr_faces/(sqrt(length(logPwrAll_faces))); %SEM of log transformed power
        
        figure(22)
        if lowIdx == 1 && PFCIdx == 1 || lowIdx == 1 && PFCIdx == 2 || lowIdx == 1 && PFCIdx == 3 || ...
                lowIdx == 1 && PFCIdx == 4
            subplot(4,2,PFCIdx)
            thisPlot_PFC(:,counterP) = mLogPwr_faces(:,2); counterP = counterP+1;
            errorbar(mLogPwr_faces(:,1),semLogPwr_faces(:,1), 'LineWidth', 2); hold on;
            ylim([-5 0.2]); xlim([1 length(mCoh_faces)]); % ylim([-3 0]); %ylabel('Coherence','FontSize',12);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),length(freqWindow)], 'XTickLabel', ...
                {ceil(freqWindow(1)),ceil(freqWindow(ceil(length(freqWindow)/2))),ceil(freqWindow(end))}); xlim([1 length(mCoh_faces)]);
            set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]);
            title(['log Pwr: ' num2str(usePFCCoord)], 'FontSize', 15);
        end
        if lowIdx == 1 && PFCIdx == 1
            subplot(4,2,5)
            thisPlot_Low(:,counterL) = mLogPwr_faces(:,2); counterL = counterL+1;
            errorbar(mLogPwr_faces(:,2),semLogPwr_faces(:,2), 'LineWidth', 2); hold on;
            set(gca, 'FontSize', 15); %xlabel('Frequency (Hz)','FontSize',15);
            ylim([-5 0.2]); xlim([1 length(mCoh_faces)]); % ylim([-3 0]); %ylabel('Coherence','FontSize',12);
            title(['log Pwr: ' num2str(useLowLvlCoord)], 'FontSize', 15); %legend({'', 'TEST', '', 'TRAIN'});
            set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]);
        elseif lowIdx == 2 && PFCIdx == 2
            subplot(4,2,6)
            thisPlot_Low(:,counterL) = mLogPwr_faces(:,2); counterL = counterL+1;
            errorbar(mLogPwr_faces(:,2),semLogPwr_faces(:,2), 'LineWidth', 2); hold on;
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            ylim([-5 0.2]);  xlim([1 length(mCoh_faces)]); % xlim([0 10]) % ylim([-3 0]); %ylabel('Coherence','FontSize',12);
            title(['log Pwr: ' num2str(useLowLvlCoord)], 'FontSize', 15); %legend({'', 'TEST', '', 'TRAIN'});
            set(gca, 'XTick', [1,round(length(freqWindow)/2),length(freqWindow)], 'XTickLabel', ...
                {ceil(freqWindow(1)),ceil(freqWindow(ceil(length(freqWindow)/2))),ceil(freqWindow(end))}); xlim([1 length(mCoh_faces)]);
        elseif lowIdx == 3 && PFCIdx == 3
            subplot(4,2,7)
            thisPlot_Low(:,counterL) = mLogPwr_faces(:,2); counterL = counterL+1;
            errorbar(mLogPwr_faces(:,2),semLogPwr_faces(:,2), 'LineWidth', 2); hold on;
            set(gca, 'FontSize', 15); %xlabel('Frequency (Hz)','FontSize',15);
            ylim([-5 0.2]);  xlim([1 length(mCoh_faces)]); % xlim([0 10]) % ylim([-3 0]); %ylabel('Coherence','FontSize',12);
            title(['log Pwr: ' num2str(useLowLvlCoord)], 'FontSize', 15); %legend({'', 'TEST', '', 'TRAIN'});
            set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]);
        elseif lowIdx == 4 && PFCIdx == 4
            subplot(4,2,8)
            thisPlot_Low(:,counterL) = mLogPwr_faces(:,2); counterL = counterL+1;
            errorbar(mLogPwr_faces(:,2),semLogPwr_faces(:,2), 'LineWidth', 2); hold on;
            set(gca, 'FontSize', 15); %xlabel('Frequency (Hz)','FontSize',15);
            ylim([-5 0.2]);  xlim([1 length(mCoh_faces)]); % xlim([0 10]) % ylim([-3 0]); %ylabel('Coherence','FontSize',12);
            title(['log Pwr: ' num2str(useLowLvlCoord)], 'FontSize', 15); %legend({'', 'TEST', '', 'TRAIN'});
            set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Load Random Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load(loadNameRandom); %load random data
        
        for ss = 1:size(allPowerInfo,1)
            pwrData_random = allPowerInfo(ss).powspctrm(:,1:length(freqWindow))
            pwrAll_random(ss,:,:) = pwrData_random'; % power (form trials x freq x pwrPerCond)
        end
        
        mPwr_random = squeeze(mean(pwrAll_random,1))
        stdPwr_random = squeeze(std(pwrAll_random))
        semPwr_random = stdPwr_faces/(sqrt(length(pwrAll_random)))
        
        %Log power
        logPwrAll_random = log10(pwrAll_random*10);
        mLogPwr_random = squeeze(mean(logPwrAll_random,1))
        stdLogPwr_random = squeeze(std(logPwrAll_random))
        semLogPwr_random = stdLogPwr_random/(sqrt(length(pwrAll_random)))
        
        figure(22)
        if lowIdx == 1 && PFCIdx == 1 || lowIdx == 1 && PFCIdx == 2 || lowIdx == 1 && PFCIdx == 3 || ...
                lowIdx == 1 && PFCIdx == 4
            subplot(4,2,PFCIdx)
            errorbar(mLogPwr_random(:,1),semLogPwr_random(:,1), 'LineWidth', 2); hold on;
        end
        if lowIdx == 1 && PFCIdx == 1
            subplot(4,2,5)
            errorbar(mLogPwr_random(:,2),semLogPwr_random(:,2), 'LineWidth', 2); hold on;
        elseif lowIdx == 2 && PFCIdx == 2
            subplot(4,2,6)
            errorbar(mLogPwr_random(:,2),semLogPwr_random(:,2), 'LineWidth', 2); hold on;
        elseif lowIdx == 3 && PFCIdx == 3
            subplot(4,2,7)
            errorbar(mLogPwr_random(:,2),semLogPwr_random(:,2), 'LineWidth', 2); hold on;
        elseif lowIdx == 4 && PFCIdx == 4
            subplot(4,2,8)
            errorbar(mLogPwr_random(:,2),semLogPwr_random(:,2), 'LineWidth', 2); hold on;
        end
        
        
        %All power data
        all_Logpwr_PFC = [mLogPwr_faces(:,1)'; mLogPwr_random(:,1)']
        all_Logpwr_PFC_sem = [semLogPwr_faces(:,1)'; semLogPwr_random(:,1)']
        all_Logpwr_FFA = [mLogPwr_faces(:,2)'; mLogPwr_random(:,2)']
        all_Logpwr_FFA_sem = [semLogPwr_faces(:,2)'; semLogPwr_random(:,2)']
        
        if lowIdx == 1 && PFCIdx == 1 || lowIdx == 2 && PFCIdx == 2 || ...
                lowIdx == 3 && PFCIdx == 3 || lowIdx == 4 && PFCIdx == 4
            figure(10022)
            subplot(2,2,countBB)
            %Plot for FFA
            mseb((repmat(1:length(all_Logpwr_FFA),2,1)),all_Logpwr_FFA,all_Logpwr_FFA_sem,[],1); box on;
            set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]); xlim([1 length(all_m_Coh)]); ylim([-2.6 0.1]);
            countBB = countBB + 1;
        end
        %Only plot PFC power in particular loops (otherwise double up
        %plots)
        if lowIdx == 1 && PFCIdx == 1 || lowIdx == 1 && PFCIdx == 2 || lowIdx == 1 && PFCIdx == 3 || ...
                lowIdx == 1 && PFCIdx == 4
            figure(10033)
            subplot(2,2,countB)
            %Plot for PFC
            mseb((repmat(1:length(all_Logpwr_PFC),2,1)),all_Logpwr_PFC,all_Logpwr_PFC_sem,[],1); box on;
            if countB == 3
                set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
                set(gca, 'XTick', [1,round(length(freqWindow)/2),length(freqWindow)], 'XTickLabel', ...
                    {ceil(freqWindow(1)),ceil(freqWindow(ceil(length(freqWindow)/2))),ceil(freqWindow(end))});
                xlim([1 length(all_m_Coh)]); ylim([-2.6 0.1]);
            else
                set(gca,'YTickLabel',[]);  set(gca,'XTickLabel',[]); xlim([1 length(all_m_Coh)]); ylim([-2.6 0.1]);
            end
        end
        countB = countB + 1;
        
        
        
        %% PLOT IN ONE FIGURE (for comparison with simulated results)
           figure(88)
           subplot(1,5,2)
            %Plot POWER for PFC
            mseb((repmat(1:length(all_Logpwr_PFC),2,1)),all_Logpwr_PFC,all_Logpwr_PFC_sem,[],1)
          set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),ceil(length(freqWindow))], 'XTickLabel', {ceil(freqWindow(1)),...
                freqWindow(ceil(length(freqWindow)/2)),freqWindow(ceil(length(freqWindow)))}); xlim([1 length(mCoh_faces)]);
            title(['Power L DPFC'], 'FontSize', 15); 
           figure(88)
           subplot(1,5,1)
         
            %Plot POWER for FFA
            mseb((repmat(1:length(all_Logpwr_FFA),2,1)),all_Logpwr_FFA,all_Logpwr_FFA_sem,[],1)
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),ceil(length(freqWindow))], 'XTickLabel', {ceil(freqWindow(1)),...
                freqWindow(ceil(length(freqWindow)/2)),freqWindow(ceil(length(freqWindow)))}); xlim([1 length(mCoh_faces)]);
            title(['Power L FG'], 'FontSize', 15); 
         
            %Plot GC to PFC
            figure(88)
           subplot(1,5,4)
            mseb((repmat(1:length(all_GC_to_PFC),2,1)),all_GC_to_PFC,all_GC_to_PFC_sem,[],1)
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),ceil(length(freqWindow))], 'XTickLabel', {ceil(freqWindow(1)),...
                freqWindow(ceil(length(freqWindow)/2)),freqWindow(ceil(length(freqWindow)))}); xlim([1 length(mCoh_faces)]);
            title(['GC L FG --> L DPFC'], 'FontSize', 15); 

            %Plot GC from PFC
                    figure(88)
           subplot(1,5,5)
            mseb((repmat(1:length(all_GC_from_PFC),2,1)),all_GC_from_PFC,all_GC_from_PFC_sem,[],1)
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),ceil(length(freqWindow))], 'XTickLabel', {ceil(freqWindow(1)),...
                freqWindow(ceil(length(freqWindow)/2)),freqWindow(ceil(length(freqWindow)))}); xlim([1 length(mCoh_faces)]);
            title(['GC L DPFC--> L FG'], 'FontSize', 15); 
         
            %Plot coherence
            figure(88)
           subplot(1,5,3)
            mseb((repmat(1:length(all_m_Coh),2,1)),all_m_Coh,all_sem_Coh,[],1)
            set(gca, 'FontSize', 15); xlabel('Frequency (Hz)','FontSize',15);
            set(gca, 'XTick', [1,round(length(freqWindow)/2),ceil(length(freqWindow))], 'XTickLabel', {ceil(freqWindow(1)),freqWindow(ceil(length(freqWindow)/2)),freqWindow(ceil(length(freqWindow)))}); xlim([1 length(mCoh_faces)]);
            title(['Coherence'], 'FontSize', 15); 
            
        
        
        
    end
end
