%% Plot power spectrum for Face and Random single-trials (mean across trials)

% This script uses the Chronux toolbox for power spectrum estimation
%   using mtspecramc

% This script builds on Steps 1 and 2 for non-parametric GC

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
% STEP 3: Plot power spectrum by subtracting baseline power at each freq.

% Elise Rowe, Monash University, 2020

clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT: Settings and filenames for extracting ROI data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filepath = ['/Users/egrow1/Desktop/ALL_SW_0_to_500ms_EBRem/']; %set file path
filepathB = ['/Users/egrow1/Desktop/ALL_SW_noF_0_to_500ms_EBR/']; %set file path

filepathList = ['/Users/egrow1/Desktop/NEW_Redo_OriginalChans_SHAFTO/SourceResults_MASSIVE/'];
filePathCoords = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/1_EEG_PreProcessing/'];
count = 0;
% SELECT THE TWO NODES TO UNDERGO GRANGER CAUSALITY ANALYSIS
%possCoords = {'RightOcc','RightFFA','RightPar','RightITG','RightPFC'} %,'LeftPar'}
possCoords = {'L_PFC'}
% possible left coordinates: 'LeftITG','LeftOcc','LeftPar','LeftFFA'
% possible right coordinates:

nNodes = 2; %number of nodes in the GC network
nTrialsFace = 203; %number of face trials
nTrialsRandom = 348; %number of random trials
coordList = []; %empty field for coordinate list
timeWindow = [0 500]; %examine this time window
timeIdx = linspace(0,500,251); % find index for time points in data
useWind = [find(timeIdx==timeWindow(1)):find(timeIdx==timeWindow(2))]; %assign idxs by time window

%%%%% INPUTS FOR FIELDTRIP %%%%
nTimePts = 251; %number of time points
sampFreq = 500; %sampling frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1: EXTRACT SOURCEWAVE PER TRIAL FROM SINGLE-TRIAL DATAFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for gg = 1:length(possCoords)
    
    listCoords = ['Coord_NoF_FINAL_' num2str(possCoords{gg}) '_by_TstatOverTime_0_to_500ms'];
    load([filePathCoords listCoords], 'sortedCoords')
    coordList =sortedCoords(1,:)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Extract waveforms for face trials
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for voxel = 1:size(coordList,1)
        thisVoxel = coordList(voxel,:);
        filename = ['NEW_Source_Waveform_0to500ms_at_-32 -52  -6_for_ALL_trials_mSae_faces_aEBfspmeeg_7MJN_0to500ms_trial_205.mat'];
   %     filename = ['Source_Waveform_at_-32 -52  -6_for_ALL_trials_fmSae_faces_aEBfspmeeg_7MJN_0to500ms_trial_205.mat'];
     %     filename = ['NEW_Source_Waveform_0to500ms_at_-8   1  69_for_ALL_trials_mSae_faces_aEBfspmeeg_7MJN_0to500ms_trial_205.mat'];
   %     filename = ['Source_Waveform_at_-8   1  69_for_ALL_trials_fmSae_faces_aEBfspmeeg_7MJN_0to500ms_trial_205.mat'];
     load ([filepathB filename])
        GM_faces(voxel,:,:) = waveform;% feature x time x trial
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Extract waveforms for RANDOM trials
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for voxel = 1:size(coordList,1)
        thisVoxel = coordList(voxel,:);
 %       filename = ['NEW_Source_Waveform_0to500ms_at_-32 -52  -6_for_ALL_trials_mSae_random_aEBfspmeeg_7MJN_0to500ms_trial_354.mat'];
      filename = ['Source_Waveform_at_-32 -52  -6_for_ALL_trials_fmSae_faces_aEBfspmeeg_7MJN_0to500ms_trial_205.mat'];

 %     filename = ['Source_Waveform_at_-32 -52  -6_for_ALL_trials_fmSae_random_aEBfspmeeg_7MJN_0to500ms_trial_354.mat'];
 %       filename = ['NEW_Source_Waveform_0to500ms_at_-8   1  69_for_ALL_trials_mSae_faces_aEBfspmeeg_7MJN_0to500ms_trial_205.mat'];
 %        filename = ['Source_Waveform_at_-8   1  69_for_ALL_trials_fmSae_random_aEBfspmeeg_7MJN_0to500ms_trial_354.mat'];
 load ([filepath filename])
        GM_random(voxel,:,:) = waveform;% feature x time x trial
    end
    
    %% Gather all data within time window (if applicable)
    GM_faces = GM_faces(:,:,useWind); %restrict to specified time window (if applicable)
    GM_random = GM_random(:,:,useWind); %restrict to specified time window (if applicable)
    
    data_faces(gg,:,:) = squeeze(GM_faces); % extract data (trials x time)
    data_random(gg,:,:) = squeeze(GM_random); % extract data (trials x time)
end

for ff = 1:length(possCoords)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STEP 2: PROCESS TRIAL DATA TO CONFORM TO 'COVARIANCE STATIONARITY' (demean, divide std, detrend)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %EXTRACT FACE DATA PFC (node 1)
    Face_Node = squeeze(data_faces(ff,:,:)); %extract data from this node
    X_mean = mean(Face_Node); X_std = std(Face_Node); %determine mean and std
    Face_Node = Face_Node-X_mean; % Demean (i.e. subtract the mean wave across trials from this single trial)
    for gg = 1:size(Face_Node,2)
        Face_Node(:,gg) = Face_Node(:,gg)/X_std(:,gg); %divide by stdev at each timepoint
    end
    Face_Node = detrend(Face_Node); %detrend (i.e. remove data trend)
    
    % EXTRACT RANDOM DATA PFC (node 1)
    Random_Node = squeeze(data_random(ff,:,:));  %extract data from this node
    X_mean = mean(Random_Node); X_std = std(Random_Node);  %determine mean and std
    Random_Node = Random_Node-X_mean; % Demean (subtract the mean wave across trials from this single trial)
    for gg = 1:size(Random_Node,2)
        Random_Node(:,gg) = Random_Node(:,gg)/X_std(:,gg);  %divide by stdev at each timepoint
    end
    Random_Node = detrend(Random_Node); %detrend (i.e. remove data trend)

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Spectrograms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for ss = 1:2
        count = count + 1;
        
        if ss == 1
            useData = Face_Node';
            nTrials = 203; %number of trials
        elseif ss == 2
            useData = Random_Node';
            nTrials = 348; %number random trials
        end
        
        %SETUP the PARAMETERS FOR PWR SPECTRUM IMAGING
        params.movingwin = [0.1 0.05]; %moving window (in secs)
        params.tapers = [3 5]; %tapers to use
        params.pad = 0; % padding
        params.Fs = 500; %sampling rate
        
        %EXTRACT PWR SPECTRUM ACROSS EACH TRIAL
        for iTrial = 1%:nTrials
            
            tmp = useData(:,iTrial);
            
            [pwr, time, freq] = mtspecgramc(tmp,params.movingwin,params);
            t2 = time;
            
            [p, t, f] = mtspectrumc(tmp, params)
            
            vt_roi = 0 < t2 & t2 < 0.5; % set time boundary
            vf_roi = 0 < freq & freq < 40; %set the freq. boundary
            
            ls = 10*log10(pwr); %log 10 of spectrum ('s' in form time x freq)
            ls_wind = ls(vt_roi,vf_roi,:); %spectrum within t x f window
            
            bs_ls = mean(ls(vt_roi,vf_roi,:),1) ; % baseline spectrum = mean each freq
            mbs_ls = repmat( bs_ls , [ length(time) ,  1 1]); %reshape mean log power each freq
            
            dB = ls_wind - mbs_ls; %
            alldB(:,:,iTrial) = dB'; %use this for SPECTS
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% Plotting the Spectrograms (mdB) %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% TAKE MEAN ACROSS TRIALS, AND PLOT POWER SPECTRUM
        fUsed = freq(vf_roi)
        figure(88);         
        subplot(1,2,count)
        mean_alldB = nanmean(alldB,3);
        maxVal = max(mean_alldB(:)); minVal = min(mean_alldB(:));
        imagesc(t2,fUsed,mean_alldB,[-6 6]), axis xy; 
       if ff == 1 && ss == 1;  set(gca,'FontSize', 20); ylabel('Frequency (Hz)'); xlabel('Time (secs)');
       else set(gca,'YTickLabel',[], 'XTickLabel', []); end
    if ss == 1; title('F'); set(gca,'FontSize', 20); elseif ss == 2; title('R'); set(gca,'FontSize', 20); end
    
     if ss == 1   
    p1 = p
     else
         p2 = p
     end
    
        clearvars -except data_faces f data_random ss count ff Face_Node Random_Node
    end
    
end

    plot(freq, p1, 'r'), hold on
plot(freq, p2, 'b'),
legend('lowpass','raw')
    