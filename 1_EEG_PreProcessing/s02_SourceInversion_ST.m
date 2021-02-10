%% Source inversion for single-trial EEG data

% This script uses the SPM12 toolbox for source inversion

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
clear all

%% INPUT: Settings and filenames for pre-processing
filepath = ['/Users/egrow1/Desktop/NEW_Redo_OriginalChans_SHAFTO/']
nTrials = 206; %assign number of trials
epochTimeWindow = [0 600]; %time window for epoching
trialType = 'faces'; %define trial type to run: either 'faces'or 'random'

filenameFaces = ['fmSae_' num2str(trialType) '_fspmeeg_7MJN_' num2str(epochTimeWindow(1)) ...
    'to' num2str(epochTimeWindow(2)) '_trial_' num2str(ii) '.mat']; % filename for FACE trials
filenameRandom = ['fmSae_' num2str(trialType) '_fspmeeg_7MJN_' num2str(epochTimeWindow(1)) ...
    'to' num2str(epochTimeWindow(2)) 'ms_trial_' num2str(ii) '.mat']; %fielname for RANDOM
 
% LOAD all required auxiliary files
load('chanFudicials_ShaftoFINAL.mat'); %load 'fudicials' for nas, lpa and rpa positions
load('badFaceTrials_7MJNMove_Markers-edf.mat'); %load which 'bad' face trials to skip
load('badRandomTrials_7MJNMove_Markers-edf.mat'); %load which 'bad' random trials to skip


%% Run source inversion process for each trial
% Steps include: head template, coregistrations, define and invert forward
% model, specify time/frequency window and whether to use taper, create
% images of the results

for ii = 1:nTrials
    
    spm('defaults', 'EEG')
    spm_jobman('initcfg');
    
    %Set up current filename here
    if strcmp(trialType, 'faces')
         if any(badTrials_Face(:) == ii); continue; end  %skip 'bad' trials
        filename = filenameFaces; %assign filename for faces
    elseif strcmp(trialType, 'random')
         if any(badTrials_Random(:) == ii); continue; end %skip 'bad' trials
        filename = filenameRandom; %assign filename for random
    end
    
    %Load each single-trial datafile
    load([filepath filename])
    
    %If path needs to be re-written, uncomment here
    % D.path = '../PathSource'; %Rename the file path in the EEG files to match the current location
    % D.data.fname = ['../PathSource/participantFile_trial_' num2str(trialNo) '.dat'];
    % save(filename, 'D')
    
    %Template, Coregister, Forward Model
    matlabbatch{1}.spm.meeg.source.headmodel.D = {[filepath filename]}; %load this file
    matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
    matlabbatch{1}.spm.meeg.source.headmodel.comment = 'Source'; %run source reconstruction
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1; %mesh template (1 = 'yes', 0 = 'no')
    matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2; % 1=coarse, 2 = normal, 3 = fine
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas'; %position of nasion
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = fudicials(1,:);
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa'; %position of left periauricular
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = fudicials(2,:);
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa'; %position of right periauricular
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = fudicials(3,:);
    matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0; % 0 = 'no, 1 = 'yes'
    matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM'; %boundary element method forward model
    matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell'; % single shell forward model
    
    %Invert forward model and specify time (and frequency) window
    matlabbatch{2}.spm.meeg.source.invert.D = {[filepath filename]};
    matlabbatch{2}.spm.meeg.source.invert.val = 1;
    matlabbatch{2}.spm.meeg.source.invert.whatconditions.all = 1; %use all conditions
    matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.invtype = 'GS'; %Multiple sparse priors (greedy search)
    matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.woi = epochTimeWindow; %time window
    matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.foi = [0 256]; %frequency window
    matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.hanning = 0; %1 = Hanning taper at start and end of trial (0 = no taper, 1 = yes)
    matlabbatch{2}.spm.meeg.source.invert.modality = {'EEG'}; %modality = EEG
    
    %Show results within time (and frequency) window and create images
    matlabbatch{3}.spm.meeg.source.results.D = {[filepath filename]};
    matlabbatch{3}.spm.meeg.source.results.val = 1;
    matlabbatch{3}.spm.meeg.source.results.woi = epochTimeWindow; % time of interest
    matlabbatch{3}.spm.meeg.source.results.foi = [0 0]; % frequency window specify
    matlabbatch{3}.spm.meeg.source.results.ctype = 'trials'; % 'evoked' 'induced' or single 'trials'
    matlabbatch{3}.spm.meeg.source.results.space = 1; % 1=MNI or Native
    matlabbatch{3}.spm.meeg.source.results.format = 'image';
    matlabbatch{3}.spm.meeg.source.results.smoothing = 12; % mm %Smoothing mm^3
    
    spm_jobman('serial',matlabbatch);
    
    clear spm_jobman
    
end