%% Extract source reconstructed waveforms from SPM files (single-trial files)

% This script uses SPM12 and modified core scripts by Karl Friston.

% Two stages make up this script: (1) Extracting and saving 'real' rounded
% coordinates for every voxel and (2) Extracting the corresponding source
% waves from each voxel

% 'trials' need to be specified as either 'faces' or 'random' depending on
% which trial types you are extracting source waveforms from

clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(['../../PData_MASSIVE/spm12'])); %add SPM to the path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT: Settings and filenames for pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filepath = ['/Users/egrow1/Desktop/NEW_Redo_OriginalChans_SHAFTO/']; %set filepath
con = 1; % Condition = 1 (always 1 for both types of trials)
trials = {'faces'}; %name of the experimental conditions
nTrials = 204; %number of trials

saveCoordName = ['CoordIdxs_' num2str(trials{1}) '_all8196Voxels']; %how to save these coordinates to later extract

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract Voxel Coordinates for each of the 8196 voxels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([filepath 'fmSae_' num2str(trials{1}) '_fspmeeg_7MJN_0to500_trial_1.mat'])
model = D.other.inv{1,1};
vert = model.mesh.tess_mni.vert; %find voxel coordinates

coordList = vert; %extract sourewaves from these coordinates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1: Save ('centered', nearest integer) coordinates for each voxel
%   (used for later when extracting waveforms for decoding)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for jj = 1:length(coordList)
    
    PST = coordList(jj,:);
    
    filename = ['fmSae_' num2str(trials{1}) '_fspmeeg_7MJN_0to500_trial_1.mat']; % set filena,e
    
    load([filepath filename]);
    
    % D - SPM data structure (extract inverted data)
    %==========================================================================
    model = D.other.inv{1,1};
    
    % get solution and spatiotemporal basis
    %--------------------------------------------------------------------------
    J      = model.inverse.J; %J{con} = 8196 x 5 (voxels x temporal modes) CONDITIONAL EXPECTATION
    T      = model.inverse.T; %51 x 5 (time x  temporal modes ) TEMPORAL PROJECTOR
    Is     = model.inverse.Is; % 1 x 8196 (voxels) INDICES OF ACTIVE DIPOLES/VOXELS
    pst    = model.inverse.pst; % 1 x 51 (time) (PERISTIMULUS TIME)
    R2     = model.inverse.R2; % single value (VARIANCE IN SUBSPACES ACCOUNTED FOR BY MODEL (%))
    F      = model.inverse.F; % single value (LOG EVIDENCE)
    Nd     = model.inverse.Nd; % single value (total number dipoles/voxels)
    VE     = model.inverse.VE; % single value (VARIANCE EXPLAINED IN SPATIAL/TEMPORAL SUBSPACE (%))
    
    % - project J onto pst
    %--------------------------------------------------------------------------
    J      = J{con}*T'; % Gives J = 8196 x 5
    
    % - determine the coordinates within the cortical mesh
    vert   = model.mesh.tess_mni.vert; %find voxel coordinates
    
    % Find response at XYZ
    %--------------------------------------------------------------------------
    [i,js] = min(sum([vert(Is,1)-PST(1), vert(Is,2)-PST(2), vert(Is,3)-PST(3)].^2,2));   %js = Find real voxel coords (closest to input X,Y,Z coords)
    XYZ   = round(vert(Is(js),:));       % true XYZ coods for voxel (given the index)
    
    clear D;
    coordIndexes{jj} = XYZ; %save the 'aligned coordinates' that SPM pulled waveform from
end
save(saveCoordName, 'coordIndexes')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2: Begin Extraction of the Source Waveforms from EACH coordinate and EVERY trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for jj = 1:length(coordList)
    
    PST = coordList(jj,:);
    
    for ii = 1:nTrials
        
        filename = ['fmSae_' num2str(trials{1}) '_fspmeeg_7MJN_0to500_trial_' num2str(ii) '.mat']; % set filena,e

        load([filepath filename]);
        
        % D.other.inv{1,1}.inverse
        
        % D - SPM data structure (extract inverted data)
        %==========================================================================
        model = D.other.inv{1,1};
        
        % get solution and spatiotemporal basis
        %--------------------------------------------------------------------------
        J      = model.inverse.J; %J{con} = 8196 x 5 (voxels x temporal modes) CONDITIONAL EXPECTATION
        T      = model.inverse.T; %51 x 5 (time x  temporal modes ) TEMPORAL PROJECTOR
        Is     = model.inverse.Is; % 1 x 8196 (voxels) INDICES OF ACTIVE DIPOLES/VOXELS
        pst    = model.inverse.pst; % 1 x 51 (time) (PERISTIMULUS TIME)
        R2     = model.inverse.R2; % single value (VARIANCE IN SUBSPACES ACCOUNTED FOR BY MODEL (%))
        F      = model.inverse.F; % single value (LOG EVIDENCE)
        Nd     = model.inverse.Nd; % single value (total number dipoles/voxels)
        VE     = model.inverse.VE; % single value (VARIANCE EXPLAINED IN SPATIAL/TEMPORAL SUBSPACE (%))
        
        % - project J onto pst
        %--------------------------------------------------------------------------
        J      = J{con}*T'; % Gives J = 8196 x 5
        
        % - determine the coordinates within the cortical mesh
        vert   = model.mesh.tess_mni.vert; %find voxel coordinates
        
        % Find response at XYZ
        %--------------------------------------------------------------------------
        [i,js] = min(sum([vert(Is,1)-PST(1), vert(Is,2)-PST(2), vert(Is,3)-PST(3)].^2,2));   %js = Find real voxel coords (closest to input X,Y,Z coords)
        [i,jt] = max(abs(J(js,:)));     %jt = maximum resp TIMEPOINT (regardless of +ve or -ve) at XYZ
        
        % % gather response over time
        Jt    = J(js,:);                     % extract waveform over time at this voxel
        Js    = J(:,jt);                     % extract max RESPONSE from ALL voxels at THIS TIME POINT (deemed as max for this voxel)
        XYZ   = round(vert(Is(js),:));       % true XYZ coods for voxel (given the index)
        Jmax  = abs(sparse(Is,1,Js,Nd,1));   % maximum response per voxel (collapsed over time)
        
        % % gather confidence intervals
        % %----------------------------------------------------------------------
        qC  = model.inverse.qC(js).*diag(model.inverse.qV)';
        ci  = 1.64*sqrt(abs(qC));
        
        % % PLOT IF YOU WISH
        % %----------------------------------------------------------------------
        % figure(2)
        % plot(pst,Jt,pst,Jt + ci,':',pst,Jt - ci,':',...
        %     [PST PST],[-1 1]*maxJ,':',...
        %     'Color',Color)
        % hold on
        
        % title({sprintf('estimated response - condition %d',con), ...
        %     sprintf('at %i, %i, %i mm',XYZ(1),XYZ(2),XYZ(3))})
        % xlabel('time  ms')
        % axis square
        % hold off
        
        % SAVE THESE RESULTS
        % %----------------------------------------------------------------------
        waveform(ii,:) = Jt; %this waveform
        waveform_std(ii,:) = ci; %this std of waveform
        timesteps(ii,:) = pst; %timesteps
        clear D;
    end
    resultsName = ['Source_Waveform_at_' num2str(XYZ) '_for_ALL_trials_' num2str(filename)]; %filename
    save([filepath resultsName], 'XYZ', 'waveform', 'waveform_std', 'timesteps'); %save these variables
    
    clear waveform; clear waveform_std; clear timesteps;
end

