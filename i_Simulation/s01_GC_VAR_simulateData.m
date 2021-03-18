%% VAR - Generate synthetic data to run Granger Causality analysis on
%Script that generates VAR model data for Granger Causlity simulation

%Node X refers to the 'first node' in the system, and node Y refer to the
%'second node' (i.e. signal is assumed to flow from X to Y)

%Used for generating data for two different types of trials: faces ('F')
%and random ('R'). Here, signal is assumed to flow from XF to YF (but not
%from XR to YR)

%Format for the VAR models:
% XRt+1 = Xs*XRt + N(0,1)
% XFt+1 = Xs*XFt + N(0,1)
% YRt+1 = Ys*YRt + N(0,1)
% YFt+1 = Ys*YFt + XsY*XFt + N(0,1)

clear all; close all;

%Set some basic parameters
rand('seed', sum(100 * clock)); % seed random nuber generator by clock
t = 301; %timepoints
nTrials = 206; %number of trials

%A matrix connectivity
poss_a = [-0.5: 0.05: 0.5]; %weight for 'a' (= x self conn X)
poss_c = [-0.5: 0.05: 0.5]; %weight for 'c' (= x self conn X in Y data)
poss_d = [-0.5: 0.05: 0.5]; %weight for 'd' (= y self conn Y)

for ss = 1:length(poss_a)
    for  ff = 1:length(poss_c)
        for vv = 1:length(poss_d)

            %Model parameters
            this_a = poss_a(ss); % use this weight for 'a' (= x self conn X)
            this_c = poss_c(ff);  % use this weight for 'c' (= x self conn X in Y data)
            this_d = poss_d(vv); % use this weight for 'd' (= y self conn Y)
            
            savename = ['All_Data_Synthetic_SourceERP_multi_a=' num2str(this_a) '_c=' num2str(this_c) '_d=' num2str(this_d) '.mat'];

            for gg = 1:nTrials
                
                XF(1) = 0; %start timepoint XF (X face trials)
                XR(1) = 0; %start timepoint XR (X random trials)
                YR(1) = 0; %start timepoint YR (Y random trials)
                YF(1) = 0; %start timepoint YF (Y face trials)
                
                %Generate random noise (guassian ~N(0,1)) at each timepoint (diff for each)
                for i=1:t
                    XF_rand(i) = randn; %noise for XF
                    XR_rand(i) = randn; %noise for XR
                    YR_rand(i) = randn; %noise for YR
                    YF_rand(i) = randn; %noise for YF
                end
                %Generate AR time series (loop through each time point)(Timepoint 2 relies
                %on previous timepoint T1)
                for i=2:t
                    XF(i) = this_a*XF(i-1) + XF_rand(i);
                    XR(i) = this_a*XR(i-1) + XR_rand(1,i);
                    YR(i) = this_d*YR(i-1) + YR_rand(i);
                    YF(i) = this_d*YF(i-1) + this_c*XF(i-1) + YF_rand(i);
                end
                
                %Simulated data for each trial
                XF_all(gg,:) = XF;
                XR_all(gg,:) = XR;
                YR_all(gg,:) = YR;
                YF_all(gg,:) = YF;
                
            end
            
            save(savename, 'XF_all', 'XR_all', 'YR_all', 'YF_all')
            
            %PLOT (if required)
            % figure(1); plot(XR,'color',[0.8500, 0.3250, 0.0980],'Linewidth', 6); axis off
            % figure(2); plot(YR,'color',[0.8500, 0.3250, 0.0980],'Linewidth', 6); axis off
            % figure(3); plot(XF,'Linewidth', 6); axis off
            % figure(4); plot(YF,'Linewidth', 6); axis off
            clearvars -except possXs possXsY possYs nTrials t ss ff vv
        end
    end
end
