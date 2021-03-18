%% Pre-Procesing of EEG data into individual single-trial datafiles

% This script uses the SPM12 toolbox for all preprocessing steps

% Batch script for all stages of EEG data pre-processing
    % STEP 1: Convert data from EDF format to SPM readable dat/mat
    % MANUAL STEP: Assigned channel types and locations, add event triggers
    % STEP 2: Pre-processing beings: Highpass filter, mark eyeblink
    %   artefacts, epoch depending on trial type, remove trials marked with
    %   eyeblinks and with noise > 100 uV. 
    % STEP 3: Separate trials into individual files and save
    % STEP 4: Robustly average data and lowpass filter

%  NOTE: Prior to running this script, file must be converted  to EDF format in EEGLAB
%  NOTE: After Step 1, the EEG channels need to be manually assigned and trial triggers added

% Elise Rowe, Monash University, 2020

%% INPUT: Settings and filenames for pre-processing
epochTimeWindow = [0 600]; %epoch around this time window (ms)
samplingFreq = 500; %sampling frequency (Hz)
trialType = 'Random'; %either 'Faces' or 'Random' (will extract one of the two trial types)
filename = ['P0291_Move_Markers-edf.mat'] %Converted data filename
filepath = ['/Users/egrow1/Desktop/Shafto_GitHub_Scripts/1_EEG_PreProcessing/']
load('chanFudicials_ShaftoFINAL.mat'); %load 'fudicials' for source re-coregistration nas, lpa and rpa positions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1: Convert the datafile to SPM readable format (from EDF to dat/mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spm_jobman('initcfg');
clear matlabbatch

filename = ['P0291_Move_Markers-edf.mat'] %Data files

%%Convert (from EDF to dat/mat file)
matlabbatch{1}.spm.meeg.convert.dataset = {filename};
matlabbatch{1}.spm.meeg.convert.mode.continuous.readall = 1;
matlabbatch{1}.spm.meeg.convert.channels{1}.all = 'all';
matlabbatch{1}.spm.meeg.convert.outfile = '';
matlabbatch{1}.spm.meeg.convert.eventpadding = 0;
matlabbatch{1}.spm.meeg.convert.blocksize = 3276800;
matlabbatch{1}.spm.meeg.convert.checkboundary = 1;
matlabbatch{1}.spm.meeg.convert.saveorigheader = 0;
matlabbatch{1}.spm.meeg.convert.inputformat = 'autodetect';

spm_jobman('run',matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MANUAL PREPARATION NOW REQUIRED!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHANNELS: Manually setup the file by adding the channels using the 'Prepare' function in GUI
% (1) Manually setup the channels using the 'channel locations .xyz file'

% Need to go into this file and convert the X Y Z coordinates to two
% separate files ---- (1) for EEG channels and (2) the other for EOG channels

% MANUALLY enter into D.channels structure if not automatically in last
% step -- need to also select 'uV' as units

%% TRIGGERS: Manually need to add the triggers
% For this you need to use an old SPM files D.trials.events structure and
% enter your own labels that correspond to the information given in the
% vmrk file! (this needs to be converted into WOrd, then Excel then into
% Matlab row/col format

% INPUT: Every trigger NAME into the 'value' column
% Run the script to replace every 'type' as "STATUS"
for ii = 5:length(D.trials.events)
    D.trials.events(ii).type = 'STATUS';
end

%INPUT: Every trigger TIME into the 'time' column
%Run the script to shift these raw time values into secs
for ii = 4:length(D.trials.events)
    currTime = D.trials.events(ii).time;
    D.trials.events(ii).time = currTime/samplingFreq;
end
% 

%% STEP 2: Data pre-processing begins here: Filter and remove artefacts
% (1) Highpass filter, (2) mark eyeblink artefacts, (3) epoch the data depending 
% on trial type (defned by trigger), (3) detect noisy  artefacts (> 100 Hz) and 
% remove trials marked with eyeblinks and noisy artefacts.

spm_jobman('initcfg');
clear matlabbatch

%Highpass filter (> 0.5 Hz)
matlabbatch{1}.spm.meeg.preproc.filter.D = {[filepath 'spmeeg_' filename]}
matlabbatch{1}.spm.meeg.preproc.filter.band = 'high';
matlabbatch{1}.spm.meeg.preproc.filter.freq = 0.5; % Highpass filter frequency
matlabbatch{1}.spm.meeg.preproc.filter.dir = 'twopass';
matlabbatch{1}.spm.meeg.preproc.filter.order = 5;
matlabbatch{1}.spm.meeg.preproc.filter.prefix = 'f'; %prefix added to the filename

%Setup HEAD MODEL for eyeblink detection (using source-level coregistration)
matlabbatch{2}.spm.meeg.source.headmodel.D = {[filepath 'fspmeeg_' filename]};
matlabbatch{2}.spm.meeg.source.headmodel.val = 1;
matlabbatch{2}.spm.meeg.source.headmodel.comment = 'Source';
matlabbatch{2}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
matlabbatch{2}.spm.meeg.source.headmodel.meshing.meshres = 2; % 1=coarse, 2 = normal, 3 = fine
matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = fudicials(1,:);
matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = fudicials(2,:);
matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = fudicials(3,:);
matlabbatch{2}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
matlabbatch{2}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
matlabbatch{2}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';

%Setup eyeblink detection (use VEOG channel to mark trials with eyeblinks)
matlabbatch{3}.spm.meeg.preproc.artefact.D = {[filepath 'fspmeeg_' filename]};
matlabbatch{3}.spm.meeg.preproc.artefact.mode = 'mark';
matlabbatch{3}.spm.meeg.preproc.artefact.badchanthresh = 0.2;
matlabbatch{3}.spm.meeg.preproc.artefact.append = true;
matlabbatch{3}.spm.meeg.preproc.artefact.methods.channels{1}.chan = 'VEOG';
matlabbatch{3}.spm.meeg.preproc.artefact.methods.fun.eyeblink.threshold = 4;
matlabbatch{3}.spm.meeg.preproc.artefact.methods.fun.eyeblink.excwin = 0;
matlabbatch{3}.spm.meeg.preproc.artefact.prefix = 'aEB';

%Epoch by trialType (either Faces or Random)
matlabbatch{4}.spm.meeg.preproc.epoch.D = {[filepath 'aEBfspmeeg_' filename]};
matlabbatch{4}.spm.meeg.preproc.epoch.trialchoice.define.timewin = epochTimeWindow;
% Select trial definitions to use (depending on 'trialType' defined at start of script)
if strcmp(trialType, 'Faces')
    matlabbatch{4}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).conditionlabel = 'Phase3_Face';
    matlabbatch{4}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventtype = 'STATUS';
    matlabbatch{4}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventvalue = ' S 31';
    matlabbatch{4}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).trlshift = 0;
    matlabbatch{4}.spm.meeg.preproc.epoch.prefix = 'e_faces_';
elseif strcmp(trialType, 'Random')
    matlabbatch{4}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).conditionlabel = 'Phase3_Random';
    matlabbatch{4}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventtype = 'STATUS';
    matlabbatch{4}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventvalue = ' S 32';
    matlabbatch{4}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).trlshift = 0;
    matlabbatch{4}.spm.meeg.preproc.epoch.prefix = 'e_random_';
end
matlabbatch{4}.spm.meeg.preproc.epoch.bc = 1;
matlabbatch{4}.spm.meeg.preproc.epoch.eventpadding = 0;


%Detect artefacts at 100uV and REMOVE trials marked with eyeblinks
if strcmp(trialType, 'Faces')
    matlabbatch{5}.spm.meeg.preproc.artefact.D = {[filepath 'e_faces_aEBfspmeeg_' filename]};
elseif strcmp(trialType, 'Random')
    matlabbatch{5}.spm.meeg.preproc.artefact.D = {[filepath 'e_random_aEBfspmeeg_' filename]};
end

matlabbatch{5}.spm.meeg.preproc.artefact.mode = 'reject'; %Reject channels according to methods below
matlabbatch{5}.spm.meeg.preproc.artefact.badchanthresh = 0.2; %bad channel threshold (default = 0.2; usually used)
matlabbatch{5}.spm.meeg.preproc.artefact.append = true;
matlabbatch{5}.spm.meeg.preproc.artefact.methods(1).channels{1}.chan = 'VEOG'; %find VEOG channels w eyeblinks
matlabbatch{5}.spm.meeg.preproc.artefact.methods(1).fun.events.whatevents.artefacts = 1; % 1 = "all"
matlabbatch{5}.spm.meeg.preproc.artefact.methods(2).channels{1}.all = 'all'; %Also rejet channels w noise above 100 uV
matlabbatch{5}.spm.meeg.preproc.artefact.methods(2).fun.threshchan.threshold = 100; %100 uV
matlabbatch{5}.spm.meeg.preproc.artefact.methods(2).fun.threshchan.excwin = 1000; %default
matlabbatch{5}.spm.meeg.preproc.artefact.prefix = 'a';

spm_jobman('run',matlabbatch);

clear spm_jobman


%% STEP 3: Extract each trial as a separate data-file
% Here, we first determine if any trials were marked as 'bad' (i.e. reject)
% using artefact detection (i.e. eyeblinks or 100uV noise). These trials
% are removed from further procesing and their labels are saved for future 
% reference.

% We then extract each 'good' trial as a separate datafile (1 per trial)

if strcmp(trialType, 'Faces')
    
    load([filepath 'ae_faces_aEBfspmeeg_' filename]); %load dataset
    nTrials_Face = length(D.trials); %total number of trials

    for bb = 1:nTrials_Face  %Find bad trials (marked to reject)
        thisTrialReject_face(bb) = D.trials(bb).bad; % ('1' = bad, '0' = good)
    end
    
    %Save this 'bad' face trials for this participant (for use later)
    badTrials_Face = find(thisTrialReject_face == 1)
    save(['badFaceTrials_' filename], 'badTrials_Face')
    
    
    for ii = 1:nTrials_Face
        
        if any(badTrials_Face(:) == ii) %do not extract 'bad' trials
            continue
        end
        
        load([filepath 'ae_faces_aEBfspmeeg_' filename]); %load dataset
        
        D_orig.trials = D.trials; %extract trial data under separate name
        D.trials = []; % delete the original trial structure (to be overwritten)
        
        for kk = ii+1:length(D_orig.trials)
            D_orig.trials(kk).label = 'Phase_OTHER_Faces'; % re-label all trials BEFORE the current trial
        end
        
        if ii > 1 % re-label all other trials AFTER the current trial
            fillUp = 1:ii-1;
            for jj = 1:max(fillUp)
                D_orig.trials(jj).label = 'Phase_OTHER_Faces';
            end
        end
        
        D.trials = D_orig.trials; %re-define the trial types in this data file
        
        savename = ['Sae_faces_aEBfspmeeg_P0291__' num2str(epochTimeWindow(1)) 'to' ...
            num2str(epochTimeWindow(2)) 'ms_trial_' num2str(ii) '.mat']
        save(savename, 'D'); % save this trial as a separate file
        
        clear D_orig.trials; clear D; % clear for next loop
    end
    
elseif strcmp(trialType, 'Random')
    % Do the same for the random trials
    load([filepath 'ae_random_aEBfspmeeg_' filename])
    nTrials_Random = length(D.trials); %nuber of trials for this trial type
    
    for bb = 1:nTrials_Random %Find bad trials (marked to reject)
        thisTrialReject_rand(bb) = D.trials(bb).bad; 
    end
    
    %Save this 'bad' face trials for this participant (for later)
    badTrials_Random = find(thisTrialReject_rand == 1)% ('1' = bad, '0' = good)
    save(['badRandomTrials_' filename], 'badTrials_Random')
    
    for ii = 1:nTrials_Random
        
        if any(badTrials_Random(:) == ii) %do not extract 'bad' trials
            continue
        end
        
        load([filepath 'ae_random_aEBfspmeeg_' filename]); %load preprocessed dataset
        
        D_orig.trials = D.trials; %extract trial data under separate name
        D.trials = []; % delete the original trial structure (to be overwritten)
        
        for kk = ii+1:length(D_orig.trials)
            D_orig.trials(kk).label = 'OTHER_Random'; % re-label all trials BEFORE the current trial
        end
        
        if ii > 1
            fillUp = 1:ii-1; %relabel all the trials AFTER the current trial
            for jj = 1:max(fillUp)
                D_orig.trials(jj).label = 'OTHER_Random';
            end
        end
        
        D.trials = D_orig.trials; %re-define the trial types in this data file
        
        savename = ['Sae_random_aEBfspmeeg_P0291__' num2str(epochTimeWindow(1)) 'to' ...
            num2str(epochTimeWindow(2)) 'ms_trial_' num2str(ii) '.mat']
        save([filepath savename], 'D'); % save this trial as a separate file
        
        clear D_orig.trials; clear D; %clear for next loop
    end
end


%% STEP 4: FINALISE the pre-processing after trials have been separated into single files
% Final steps here are to robustly average (= 2 final trial types, the one
% of interest and the mean of all other trials, not used). As robust
% averaging can introduce high frequency noise, we then low-pass filter <  40 Hz.

%Run for face trials
if strcmp(trialType, 'Faces')
    for ii = 1:nTrials_Face
        
        if any(badTrials_Face(:) == ii) % skip 'bad' trials
            continue
        end
        
        savename = ['Sae_faces_aEBfspmeeg_P0291__' num2str(epochTimeWindow(1)) 'to' ...
            num2str(epochTimeWindow(2)) 'ms_trial_' num2str(ii) '.mat']
        
        spm_jobman('initcfg');
        
        %Robustly average (all defaults -- usually don't change much here)
        matlabbatch{1}.spm.meeg.averaging.average.D(1) = {[filepath savename]};
        matlabbatch{1}.spm.meeg.averaging.average.userobust.robust.ks = 3;
        matlabbatch{1}.spm.meeg.averaging.average.userobust.robust.bycondition = true;
        matlabbatch{1}.spm.meeg.averaging.average.userobust.robust.savew = false;
        matlabbatch{1}.spm.meeg.averaging.average.userobust.robust.removebad = false;
        matlabbatch{1}.spm.meeg.averaging.average.plv = false;
        matlabbatch{1}.spm.meeg.averaging.average.prefix = 'm';
        
        %Lowpass filter (robust averaging introduces high frequency noise)
        matlabbatch{2}.spm.meeg.preproc.filter.D(1) = {[filepath 'm' savename]};
        matlabbatch{2}.spm.meeg.preproc.filter.type = 'butterworth';
        matlabbatch{2}.spm.meeg.preproc.filter.band = 'low';
        matlabbatch{2}.spm.meeg.preproc.filter.freq = 40; % 40 Hz
        matlabbatch{2}.spm.meeg.preproc.filter.dir = 'twopass'; %default
        matlabbatch{2}.spm.meeg.preproc.filter.order = 5;
        matlabbatch{2}.spm.meeg.preproc.filter.prefix = 'f';
        
        spm_jobman('run',matlabbatch);
        clear spm_jobman
    end
    
% Run for random trials    
elseif strcmp(trialType, 'Random')
    for ii = 1:nTrials_Random
        
        if any(badTrials_Random(:) == ii)
            continue
        end
        
        spm_jobman('initcfg');
        
        savename = ['Sae_random_aEBfspmeeg_P0291__' num2str(epochTimeWindow(1)) 'to' ...
            num2str(epochTimeWindow(2)) 'ms_trial_' num2str(ii) '.mat']
        
        %Robustly average (all defaults -- usually don't change much here)
        matlabbatch{1}.spm.meeg.averaging.average.D(1) = {[filepath savename]};
        matlabbatch{1}.spm.meeg.averaging.average.userobust.robust.ks = 3;
        matlabbatch{1}.spm.meeg.averaging.average.userobust.robust.bycondition = true;
        matlabbatch{1}.spm.meeg.averaging.average.userobust.robust.savew = false
        matlabbatch{1}.spm.meeg.averaging.average.userobust.robust.removebad = false;
        matlabbatch{1}.spm.meeg.averaging.average.plv = false;
        matlabbatch{1}.spm.meeg.averaging.average.prefix = 'm';
        
        %Lowpass filter (robust averaging introduces high frequency noise)
        matlabbatch{2}.spm.meeg.preproc.filter.D(1) = {[filepath 'm' savename]};
        matlabbatch{2}.spm.meeg.preproc.filter.type = 'butterworth';
        matlabbatch{2}.spm.meeg.preproc.filter.band = 'low';
        matlabbatch{2}.spm.meeg.preproc.filter.freq = 40; % 40 Hz
        matlabbatch{2}.spm.meeg.preproc.filter.dir = 'twopass'; %default
        matlabbatch{2}.spm.meeg.preproc.filter.order = 5;
        matlabbatch{2}.spm.meeg.preproc.filter.prefix = 'f';
        
        spm_jobman('run',matlabbatch);
        clear spm_jobman
        
    end
end