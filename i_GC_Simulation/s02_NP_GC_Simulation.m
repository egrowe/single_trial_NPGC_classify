%% Extract individual ERPs and plot over time
clear all; close all;
%Filepath for loading datafiles
filepath = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/i_GC_Simulation/'];
%Add fieldtrip and bsmart toolbox to the path
addpath(genpath('/Users/egrow1/Desktop/Shafto_Analysis_Scripts/Single-Trial/Decoding/bsmart'));
addpath('/Users/egrow1/Desktop/spm12/external/fieldtrip'); ft_defaults;

%SETUP THE PARAMETERS
%A matrix connectivity
poss_a = [0.5]; %weight for 'a' (= x self conn X)
poss_c = [-0.5, 0, 0.5]; %weight for 'c' (= x self conn X in Y data)
poss_d = [-0.5, 0, 0.5];  %weight for 'd' (= y self conn Y)

%Setup some other major parameters
theseTrials = 'faces'; %options = 'random' or 'faces': analyse particular trial data
nNodes = 2; %number of nodes in GC calculation
nTrials = 206; %total number of trials (note, each loded individually for single-trial analysis)

%Assign other parameters
timeWindow = [0 500]; %examine this time window
nTimePts = 251; %number of timepoints within each trial
sampRate = 500; %sampling rate (Hz)

count = 0;

for gg = 1:length(poss_a)
    
    for  ff = 1:length(poss_c)
        
        for vv = 1:length(poss_d)
            count = count + 1;

            %Model parameters
            this_a = poss_a(gg); % use this weight for 'a' (= x self conn X)
            this_c = poss_c(ff);  % use this weight for 'c' (= x self conn X in Y data)
            this_d = poss_d(vv); % use this weight for 'd' (= y self conn Y)

            %Load this data set
            loadname = ['All_Data_Synthetic_SourceERP_multi_Xs=' num2str(this_a) '_XsY=' num2str(this_c) '_Ys=' num2str(this_d) '.mat'];
            load([loadname])
            
            %Save name
            savename = ['all' num2str(theseTrials) '_CohPower_GC_SYNTHETIC_XX=' num2str(this_a) '_XY=' num2str(this_c) '_YY=' num2str(this_d) '.mat']
 
            %EXTRACT DATA
            Face_PFC = YF_all(:,1:nTimePts); % Face data at PFC (or 'first node', i.e 'YF')
            Face_FG = XF_all(:,1:nTimePts); %Face data at FG (or 'second node' i.e. 'XF')
            Random_PFC = YR_all(:,1:nTimePts); %Random data at PFC (or 'first node', i.e. 'YR')
            Random_FG =XR_all(:,1:nTimePts); %Random data at FG (or 'second' node, i.e. 'XR')
            
            %% Run GC 2 node analysis (calcuated on single-trial basis)
            for ii = 1:nTrials
                
                % Put data into fieldtrip format
                time = linspace(0,timeWindow(2),nTimePts)/1000; %time point indices
                
                for gg = 1
                    data.time{gg} = time
                end
                
                if strcmp(theseTrials, 'random') %extract trial 1, 2, 3 one-by-one each loop (if 'random')
                    for dd = 1
                        data.trial{1,dd}(1,:) = Random_FG(ii,:); %trial X data for first node
                        data.trial{1,dd}(2,:) = Random_PFC(ii,:); %trial X data for second node
                    end
                    data.label{1,:} = 'Random_Other'; %trial labels
                    data.label{2,:} = 'Random_PFC'; %trial labels
                elseif strcmp(theseTrials, 'faces') %extract trial 1, 2, 3 one-by-one each loop (if 'faces')
                    for dd = 1
                        data.trial{1,dd}(1,:) = Face_FG(ii,:); %trial X data for first node
                        data.trial{1,dd}(2,:) = Face_PFC(ii,:); %trial X data for second node
                    end
                    data.label{1,:} = 'Face_Other'; %trial labels
                    data.label{2,:} = 'Face_PFC'; %trial labels
                end
                
                data.fsample = sampRate; % sampling freq
                
                %% Computation of the multivariate autoregressive model (parametric)
                cfg         = [];
                cfg.ntrials     = nTrials; %number of trials
                cfg.triallength = 1; %length of trial
                cfg.fsample     = data.fsample; %sampling frequency
                cfg.nsignal     = 2; %number of signals
                cfg.method      = 'ar'; %method = auto-regressive
                cfg.order   = 1; %model order = 1
                cfg.method = 'bsmart'; %use toolbox
                mdata       = ft_mvaranalysis(cfg, data); %apply multivariate analysis
                
                cfg.method = 'mvar';
                mfreq      = ft_freqanalysis(cfg, mdata); %obtain spectral estimates
                
                %% NON-PARAMETRIC POWER CALCULATION
                cfg.method    = 'mtmfft'; %multitaper fourier transform
                cfg.taper     = 'dpss'; %taper type = multi-taper DPSS
                cfg.output    = 'fourier'; %output in fourier transform
                cfg.tapsmofrq = 3; %number of tapers (i.e. smoothing window, here, +/-3 Hz = 6 Hz)
                freq_NP          = ft_freqanalysis(cfg, data); %apply nonparetric frequency analysis
                fd_NP            = ft_freqdescriptives(cfg, freq_NP); %gather power descriptives
                
                %% Computation and inspection of the connectivity measures
                cfg.method    = 'coh';
                coherence_NP     = ft_connectivityanalysis(cfg, freq_NP); %calculate coherence measures
                
                %% Computation of GRANGER
                cfg           = [];
                cfg.method    = 'granger'; %apply Granger causality measures
                granger       = ft_connectivityanalysis(cfg, freq_NP);                
                
                %% GATHER ALL DATA
                allPowerInfo(ii,:) = fd_NP; %gather power data
                allCoherenceInfo(ii,:) = coherence_NP; %gather coherence data
                allGranger(ii,:) = granger; %gather GC data
                mFreq(ii,:) = mfreq; %gather data relevant to transfer/spectral density
                
                close all
                clear fd_NP; clear coherence_NP; clear granger; clear mfreq
 
            end
            
            save(savename, 'allPowerInfo','allCoherenceInfo', 'allGranger', 'mFreq')
            clear allPowerInfo; clear allCoherenceInfo; clear allGranger; clear mFreq;

        end
    end
end