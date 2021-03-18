%% Non-Parametric Granger Causality analysis on single-trial (trial-by-trial basis)

% This script requires the FieldTrip and bsmart toolbox for Granger Causality estimation

% This script is adapted from those supplied on the FieldTrip Connectivity
% Tutorial (http://www.fieldtriptoolbox.org/tutorial/connectivity/)

% Batch script for all stages of single-trial (i.e. trial-by-trial estimated)
% Non-Parametric Granger Causality estimation

% INPUT: Assign all input variables here
% STEP 1: Extract the source waveform data from the 'top voxel' (i.e.
%   best voxel according to t-stat comparison of face v random trials over time)
% STEP 2: Process trial data (individual trial level) to conform to
%   'covariance stationarity' that is required for Granger Causality
%    input data. Steps here include (1) 'demeaning' = subtracting the
%    mean across all trials from the individual trial, (2) dividing by
%    the standard deviation at each time point (using the std dev across
%    all trials at that time point) and (3) 'detrending' = remove the
%    overall trend of the data from each individual single-trial
% STEP 3: Apply non-parametric Granger Causality analysis on a
%   trial-by-trial (single-trial) basis. Outputs from this step include
%   (1) 'allGranger' = Granger coefficients from regression analysis,
%   (2) 'allPower' = Signal power at that node (at each frequency band)
%   (3) 'allCoherence' = Phase coherence between the two nodes (at each freq. band)

% If you wish to visualise some of the Granger Causality results;
%   you can uncomment the plotting code throughout Step 3.

% Elise Rowe, Monash University, 2020

clear all; close all;
addpath(genpath('/Users/egrow1/Desktop/Shafto_Analysis_Scripts/Single-Trial/Decoding/bsmart'));
addpath(genpath('/Users/egrow1/Desktop/spm12'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT: Settings and filenames for extracting ROI data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filepath = ['/Users/egrow1/Desktop/ALL_SourceWaves_8196_EBRemoved/']; %set file path
filepathList = ['/Users/egrow1/Desktop/NEW_Redo_OriginalChans_SHAFTO/SourceResults_MASSIVE/'];
filePathCoords = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/1_EEG_PreProcessing/'];

useCoordPairs = {'L_PFC','L_Occ'; 'L_PFC','L_FFA'; ...
    'R_Inf_PFC','R_Par'; 'R_Inf_PFC','R_Occ';'R_Inf_PFC','R_FFA'; ...
    'R_PFC','R_Par';'R_PFC','R_Occ'; 'R_PFC','R_FFA';...
    'L_Inf_PFC','R_Occ';'L_Inf_PFC','R_FFA'; ...
    'L_PFC','R_Par';'L_PFC','R_Occ'; 'L_PFC','R_FFA'; ...
    'R_Inf_PFC','L_Par'; 'R_Inf_PFC','L_Occ';'R_Inf_PFC','L_FFA'; ...
    'R_PFC','L_Par';'R_PFC','L_Occ'; 'R_PFC','L_FFA'}

for gg = 1:size(useCoordPairs,1)
    % SELECT THE TWO NODES TO UNDERGO GRANGER CAUSALITY ANALYSIS
    possCoords = {useCoordPairs{gg,1},useCoordPairs{gg,2}}
    
    nNodes = 2; %number of nodes in the GC network
    theseTrials = 'random'; %'random' or 'faces'
    nTrialsFace = 203; %number of face trials
    nTrialsRandom = 348; %number of random trials
    coordList = []; %empty field for coordinate list
    timeWindow = [0 600]; %examine this time window
    timeIdx = linspace(0,600,301); % find index for time points in data
    useWind = [find(timeIdx==timeWindow(1)):find(timeIdx==timeWindow(2))]; %assign idxs by time window
    
    savename = ['GCest_' num2str(theseTrials) '_' num2str(possCoords{1}) ...
        '_and_' num2str(possCoords{2}) '_0to500ms_4Tapsmo_TopVoxel.mat']; %savename
    
    %%%%% INPUTS FOR FIELDTRIP %%%%
    nTimePts = 301; %number of time points
    sampFreq = 500; %sampling frequency
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Set nTrials: Depending on the type of trials we are analysing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(theseTrials, 'faces')
        nTrials = nTrialsFace
    elseif strcmp(theseTrials, 'random')
        nTrials = nTrialsRandom
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP 1: EXTRACT SOURCEWAVE PER TRIAL FROM SINGLE-TRIAL DATAFILES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for hh = 1:length(possCoords)
        
        listCoords = ['CoordSortedEBRem_' num2str(possCoords{hh}) '_by_TstatOverTime_0_to_500ms'];
        load([filePathCoords listCoords], 'sortedCoords')
        coordList =sortedCoords(1,:)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Extract waveforms for face trials
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for voxel = 1:size(coordList,1)
            thisVoxel = coordList(voxel,:);
            filename = ['SourceWaveform_at_' num2str(thisVoxel) '_allTrials_mSae_faces_aEBfspmeeg_P0291_0to500ms.mat'];
            load ([filepath filename])
            GM_faces(voxel,:,:) = waveform;% feature x time x trial
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Extract waveforms for RANDOM trials
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for voxel = 1:size(coordList,1)
            thisVoxel = coordList(voxel,:);
            filename = ['SourceWaveform_at_' num2str(thisVoxel) '_allTrials_mSae_random_aEBfspmeeg_P0291_0to500ms.mat'];
            load ([filepath filename])
            GM_random(voxel,:,:) = waveform;% feature x time x trial
        end
        
        %% Gather all data within time window (if applicable)
        GM_faces = GM_faces(:,:,useWind); %restrict to specified time window (if applicable)
        GM_random = GM_random(:,:,useWind); %restrict to specified time window (if applicable)
        
        data_faces(hh,:,:) = squeeze(GM_faces); % extract data (trials x time)
        data_random(hh,:,:) = squeeze(GM_random); % extract data (trials x time)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP 2: PROCESS TRIAL DATA TO CONFORM TO 'COVARIANCE STATIONARITY' (demean, divide std, detrend)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %EXTRACT FACE DATA PFC (node 1)
    Face_PFC = squeeze(data_faces(1,:,:)); %extract data from this node
    X_mean = mean(Face_PFC); X_std = std(Face_PFC); %determine mean and std
    Face_PFC = Face_PFC-X_mean; % Demean (i.e. subtract the mean wave across trials from this single trial)
    for gg = 1:size(Face_PFC,2)
        Face_PFC(:,gg) = Face_PFC(:,gg)/X_std(:,gg); %divide by stdev at each timepoint
    end
    Face_PFC = detrend(Face_PFC); %detrend (i.e. remove data trend)
    
    % EXTRACT FACE DATA FFA (node 2)
    Face_Other =squeeze(data_faces(2,:,:));  %extract data from this node
    XR_mean = mean(Face_Other); XR_std = std(Face_Other);  % determine mean and std
    Face_Other = Face_Other-XR_mean; % Demean (subtract the mean wave across trials from this single trial)
    for gg = 1:size(Face_Other,2)
        Face_Other(:,gg) = Face_Other(:,gg)/XR_std(:,gg); %divide by stdev at each timepoint
    end
    Face_Other = detrend(Face_Other); %detrend (i.e. remove data trend)
    
    % EXTRACT RANDOM DATA PFC (node 1)
    Random_PFC = squeeze(data_random(1,:,:));  %extract data from this node
    X_mean = mean(Random_PFC); X_std = std(Random_PFC);  %determine mean and std
    Random_PFC = Random_PFC-X_mean; % Demean (subtract the mean wave across trials from this single trial)
    for gg = 1:size(Random_PFC,2)
        Random_PFC(:,gg) = Random_PFC(:,gg)/X_std(:,gg);  %divide by stdev at each timepoint
    end
    Random_PFC = detrend(Random_PFC); %detrend (i.e. remove data trend)
    
    % EXTRACT RANDOM DATA FFA (node 2)
    Random_Other =squeeze(data_random(2,:,:));  %extract data from this node
    XR_mean = mean(Random_Other); XR_std = std(Random_Other); %divide by stdev at each timepoint
    Random_Other = Random_Other-XR_mean;  % Demean (i.e. subtract the mean wave across trials from this single trial)
    for gg = 1:size(Random_Other,2)
        Random_Other(:,gg) = Random_Other(:,gg)/XR_std(:,gg);  %divide by stdev at each timepoint
    end
    Random_Other = detrend(Random_Other); %detrend (i.e. remove data trend)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run Non-Parametric Granger Causality Analysis on a trial-by-trial basis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii = 1%:nTrials
        %% Put data into fieldtrip format
        data.time{1} = linspace(0,500,nTimePts)/1000; %time point indices
        data.fsample = sampFreq; % sampling freq
        
        if strcmp(theseTrials, 'random') %if 'random' selected, run on random trials
            data.trial{1,1}(1,:) = Random_Other(ii,:); %data node 1 (i.e. lower in hierarchy)
            data.trial{1,1}(2,:) = Random_PFC(ii,:); % data node 2 (i.e. higher in hierarchy)
            data.label{1,:} = 'Random_Other'; %label for node 1
            data.label{2,:} = 'Random_PFC'; %label for node 2
            
        elseif strcmp(theseTrials, 'faces') %if 'faces' selected, run on face trials
            data.trial{1,1}(1,:) = Face_Other(ii,:)
            data.trial{1,1}(2,:) = Face_PFC(ii,:)
            data.label{1,:} = 'Face_Other'
            data.label{2,:} = 'Face_PFC'
            
            %% Computation of the multivariate autoregressive model (parametric)
            cfg         = [];
            cfg.ntrials     = nTrials; %number of trials
            cfg.triallength = 1; %length of trial
            cfg.fsample     = data.fsample; %sampling frequency
            cfg.nsignal     = 2; %number of signals
            cfg.method      = 'ar'; %method = auto-regressive
            cfg.order   = 1; %model order = 1
            cfg.toolbox = 'bsmart'; %use toolbox
            mdata       = ft_mvaranalysis(cfg, data); %apply multivariate analysis
            
            cfg.method = 'mvar';
            mfreq      = ft_freqanalysis(cfg, mdata); %obtain spectral estimates
            
            %% NON-PARAMETRIC POWER CALCULATION
            cfg.method    = 'mtmfft'; %multitaper fourier transform
            cfg.taper     = 'dpss'; %taper type = multi-taper DPSS
            cfg.output    = 'fourier'; %output in fourier transform
            cfg.tapsmofrq = 5; %number of tapers (i.e. smoothing window)
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
            
            close all
            clear Random_PFC Random_Other Face_PFC Face_Other fd_NP coherence_NP granger
        end
        
        save(savename, 'allCoherenceInfo','allPowerInfo', 'allGranger')
        
        clear allCoherenceInfo; clear allPowerInfo; clear allGranger;
    end