%% Extract individual ERPs and plot over time
clear all; close all;
addpath(genpath('/Users/egrow1/Desktop/NEW_Redo_OriginalChans_SHAFTO/decode_libSVM')); %add decoder to the path list
addpath(genpath('/Users/egrow1/Desktop/Shafto_Analysis_Scripts/Single-Trial/Decoding/bsmart'));
addpath(genpath('/Users/egrow1/Desktop/spm12'));

%filepath = ['/Users/egrow1/Desktop/Shafto_Analysis_Scripts/Single-Trial/Decoding/']; %data filepath
filepath = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/i_GC_Simulation/'];

possXX = 0.5 %[0.1, 0.25, 0.5]
possXY = 0.5 %[0.1, 0.25, 0.5]
possYY = 0.5 %[0.1, 0.25, 0.5]
count = 0;

%Assign other parameters
timeWindow = [0 500]; %examine this time window
nTimePts = 251; %number of timepoints within each trial
sampRate = 500; %sampling rate (Hz)


for gg = 1:length(possXX)
    
    for  ff = 1:length(possXY)
        
        for vv = 1:length(possYY)
            count = count + 1;
            %Model parameters
            Xs = possXX(gg); %x self conn X
            XsY = possXY(ff); %x self conn X in Y data
            Ys = possYY(vv); %y self conn Y
            
           % loadname = ['DATA_Synthetic_3TrialBins_Faces_Random_XF_XR_YF_YR.mat'];
            loadname = ['All_Data_Synthetic_SourceERP_multi_Xs=' num2str(Xs) '_XsY=' num2str(XsY) '_Ys=' num2str(Ys) '.mat']
            load([filepath loadname])
            
            %Setup some major parameters
            theseTrials = 'random'; %options = 'random' or 'faces': analyse particular trial data
            nNodes = 2; %number of nodes in GC calculation
            nTrials = 206; %total number of trials (note, each loded individually for single-trial analysis)
            savename = ['all_forQHplot_' num2str(theseTrials) '_CohPower_GC_SYNTHETIC.mat']
       
            %EXTRACT DATA
            Face_PFC = YF_all(:,1:nTimePts); % Face data at PFC (or 'first node', i.e 'YF')
            Face_Other = XF_all(:,1:nTimePts); %Face data at FG (or 'second node' i.e. 'XF')
            Random_PFC = YR_all(:,1:nTimePts); %Random data at PFC (or 'first node', i.e. 'YR')
            Random_Other = XR_all(:,1:nTimePts); %Random data at FG (or 'second' node, i.e. 'XR')
         
            %Plot the data if you wish
            % figure(1); plot(Face_PFC)
            % figure(2); plot(Face_Other)
            % figure(3); plot(Random_PFC)
            % figure(4); plot(Random_Other)
            
            %% Run GC 2 node analysis (calcuated on single-trial basis)
            for ii = 1:nTrials
                
                % Put data into fieldtrip format
                time = linspace(0,timeWindow(2),nTimePts)/1000; %time point indices
                
                for gg = 1
                    data.time{gg} = time
                end
                
                if strcmp(theseTrials, 'random') %extract trial 1, 2, 3 one-by-one each loop (if 'random')
                    for dd = 1
                        data.trial{1,dd}(1,:) = Random_Other(ii,:); %trial X data for first node
                        data.trial{1,dd}(2,:) = Random_PFC(ii,:); %trial X data for second node
                    end
                    data.label{1,:} = 'Random_Other'; %trial labels
                    data.label{2,:} = 'Random_PFC'; %trial labels
                elseif strcmp(theseTrials, 'faces') %extract trial 1, 2, 3 one-by-one each loop (if 'faces')
                    for dd = 1
                        data.trial{1,dd}(1,:) = Face_Other(ii,:); %trial X data for first node
                        data.trial{1,dd}(2,:) = Face_PFC(ii,:); %trial X data for second node
                    end
                    data.label{1,:} = 'Face_Other'; %trial labels
                    data.label{2,:} = 'Face_PFC'; %trial labels
                end
                
                data.fsample = sampRate; % sampling freq
                
                %% Computation of the multivariate autoregressive model (parametric)
                cfg         = [];
                cfg.ntrials     = nTrials;
                cfg.triallength = 1;
                cfg.fsample     = data.fsample;
                cfg.nsignal     = 2;
                cfg.method      = 'ar';
                cfg.order   = 1;
                cfg.method = 'bsmart';
                mdata       = ft_mvaranalysis(cfg, data);
                
                cfg.method = 'mvar';
                mfreq      = ft_freqanalysis(cfg, mdata);
                
                %% NON-PARAMETRIC POWER CALCULATION
                cfg.method    = 'mtmfft';
                cfg.taper     = 'dpss';
                cfg.output    = 'fourier';
                cfg.tapsmofrq = 3;
                freq_NP          = ft_freqanalysis(cfg, data);
                fd_NP            = ft_freqdescriptives(cfg, freq_NP);
                
                %% Computation and inspection of the connectivity measures
                cfg.method    = 'coh';
                coherence_NP     = ft_connectivityanalysis(cfg, freq_NP);
                coherenceMvar = ft_connectivityanalysis(cfg, mfreq);
                
                %% Computation of GRANGER
                cfg           = [];
                cfg.method    = 'granger';
                granger       = ft_connectivityanalysis(cfg, freq_NP);
                
                cfg           = [];
                cfg.parameter = 'grangerspctrm';
                cfg.zlim      = [0 0.1];
                ft_connectivityplot(cfg, granger);
                
                
                %% GATHER ALL DATA
                allPowerInfo(ii,:) = fd_NP;
                allCoherenceInfo(ii,:) = coherence_NP;
                allGranger(ii,:) = granger;
                mFreq(ii,:) = mfreq;
                
                close all
                clearvars -except filepath timeIdx hh gg savename data_faces data_random ...
                    allPowerInfo allCoherenceInfo allGranger ii nTrials theseTrials timeIdx ...
                    Face_PFC Face_Other gg ff vv Random_PFC Random_Other nTimePts sampRate ...
                    possXs possXsY possYscount mFreq timeWindow nTimePts sampRate
 
            end
           save(savename, 'allPowerInfo','allCoherenceInfo', 'allGranger', 'mFreq')
           % clear allGranger; clear allPowerInfo; clear allCoherenceInfo;

        end
    end
end


%% Uncomment to plot the transfer function and spectral density matrix
maxFreq = 55;       
freq_roi =  0 <= allGranger(1).freq & allGranger(1).freq < maxFreq
freqWindow = allGranger(1).freq(freq_roi)

%Extracting transfer function matrix (per trial and freq)
for ii = 1:nTrials
    transferMat = (mFreq(ii).transfer(:,:,1:length(freqWindow)));
    H_tMat(ii,:,:,:) = transferMat.*conj(transferMat);
    std_H_tMat(ii,:,:,:) = std(H_tMat(ii,:,:,:));
end

%Gather H values
mH_tMat = squeeze(mean(H_tMat,1));
std_mH_tMat = squeeze(std(H_tMat,1));
sem_mQ_sMat = std_mH_tMat/sqrt(length(std_mH_tMat))

%Gathe log H value
log_HtMat = log10(H_tMat);
log_mH_tMat = squeeze(mean(log_HtMat,1));
log_std_mH_tMat = squeeze(std(log_HtMat,1));
log_sem_mQ_sMat = log_std_mH_tMat/sqrt(length(log_std_mH_tMat))

cc = 1;
for gg = 1:2
    for  ff = 1:2
        log_hMat(cc,:) = squeeze(log_mH_tMat(gg,ff,:))
        log_sem_hMat(cc,:) = squeeze(log_sem_mQ_sMat(gg,ff,:))
        cc = cc + 1;
    end
end

%PLOT
figure(1)
subplot(1,2,1)
mseb((repmat(1:length(log_hMat),4,1)),log_hMat,log_sem_hMat,[],1)
xlim([1 length(freqWindow)]); xlabel('Frequency'); title('log H(w)H*(w)');
set(gca,'FontSize',30); legend({'XX','YX','XY','YY'}); hold off;

% %PLOT
figure(2)
subplot(1,2,1)
errorbar(squeeze(log_mH_tMat(1,1,:)),squeeze(log_sem_mQ_sMat(1,1,:)), 'LineWidth',2); hold on;
errorbar(squeeze(log_mH_tMat(1,2,:)),squeeze(log_sem_mQ_sMat(1,2,:)), 'LineWidth',2); hold on;
errorbar(squeeze(log_mH_tMat(2,1,:)),squeeze(log_sem_mQ_sMat(2,1,:)), 'LineWidth',2); hold on;
errorbar(squeeze(log_mH_tMat(2,2,:)),squeeze(log_sem_mQ_sMat(2,2,:)), 'LineWidth',2); hold on;
xlim([1 length(freqWindow)]); xlabel('Frequency'); title('log H(w)H*(w)');
set(gca,'FontSize',30); legend({'XX','YX','XY','YY'}); hold off;


%Extract spectral density matrix (per trial and freq)
for gg = 1:nTrials
    spectMatrix = (mFreq(gg).crsspctrm(:,:,1:length(freqWindow)));
    Q_sMat(gg,:,:,:) = spectMatrix.*conj(spectMatrix);
    abs_Q_sMat(gg,:,:,:) = abs(spectMatrix); %% ADD THIS BELOW IF DECIDE TO USE IT!!!!!!!
end

%Gather data for S(w) matrix
mQ_sMat = squeeze(mean(Q_sMat,1));
std_mQ_sMat = squeeze(std(Q_sMat,1));
sem_mQ_sMat = std_mQ_sMat/sqrt(length(std_mQ_sMat))

%Log S(w) data
log_Q_sMat = log10(Q_sMat);
log_mQ_sMat = squeeze(mean(log_Q_sMat,1));
log_std_mQ_sMat = squeeze(std(log_Q_sMat,1));
log_sem_mQ_sMat = log_std_mQ_sMat/sqrt(length(log_std_mQ_sMat))

%Gather for plotting
cc = 1;
for gg = 1:2
    for  ff = 1:2
        log_sMat(cc,:) = squeeze(log_mQ_sMat(gg,ff,:))
        log_sem_sMat(cc,:) = squeeze(log_sem_mQ_sMat(gg,ff,:))
        cc = cc + 1;
    end
end

%PLOT
figure(1)
subplot(1,2,2)
mseb((repmat(1:length(log_sMat),4,1)),log_sMat,log_sem_sMat,[],1)
xlim([1 length(freqWindow)]); xlabel('Frequency'); title('log Q(w)Q*(w)');
set(gca,'FontSize',30); legend({'XX','YX','XY','YY'}); hold off;


%PLOT
figure(2)
subplot(1,2,2)
errorbar(squeeze(log_mQ_sMat(1,1,:)),squeeze(log_sem_mQ_sMat(1,1,:)), 'LineWidth',2); hold on;
errorbar(squeeze(log_mQ_sMat(1,2,:)),squeeze(log_sem_mQ_sMat(1,2,:)), 'LineWidth',2); hold on;
errorbar(squeeze(log_mQ_sMat(2,1,:)),squeeze(log_sem_mQ_sMat(2,1,:)), 'LineWidth',2); hold on;
errorbar(squeeze(log_mQ_sMat(2,2,:)),squeeze(log_sem_mQ_sMat(2,2,:)), 'LineWidth',2); hold on;
xlim([1 length(freqWindow)]);xlabel('Frequency'); title('log Q(w)Q*(w)');
set(gca,'FontSize',30); legend({'XX','YX','XY','YY'}); hold off;
