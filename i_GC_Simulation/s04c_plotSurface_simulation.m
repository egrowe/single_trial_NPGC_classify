%% Plot simulation classification accuracy using SURFACE plot
clear all; close all;
filepath = ['/Users/egrow1/Desktop/GC_Sim_MASSIVE/GC_sim_MASSIVE_Results/'];
allPoss_d = [-0.5: 0.25: 0.5]; %set d outside (to allow this to be fixed at specific value)

for kk = 1:length(allPoss_d)
    
    poss_a = [-0.5: 0.05: 0.5]; %weight for 'a' (= x self conn X)
    poss_c = [-0.5: 0.05: 0.5]; %weight for 'c' (= x self conn X in Y data)
    poss_d = allPoss_d(kk); %weight for 'd' (= y self conn Y)
    
    count = 1;
    
    for gg = 1:length(poss_a)
        
        for  ff = 1:length(poss_c)
            
            for vv = 1:length(poss_d)
                
                %Model parameters
                this_a = poss_a(gg); % use this weight for 'a' (= x self conn X)
                this_c = poss_c(ff);  % use this weight for 'c' (= x self conn X in Y data)
                this_d = poss_d(vv); % use this weight for 'd' (= y self conn Y)
                
                modelParams(count,:) = [this_a, this_c, this_d];
                
                loadname = ['Results_Decode_FvR_NonPara_GC_Xs=' num2str(this_a) '_XsY=' num2str(this_c) '_Ys=' num2str(this_d) '.mat'];
                load([filepath loadname])
                
                results_NP = results; clear results; clear data;
                
                %Extract the GC decoding results to plot
                allGC_toPFC(count) = results_NP.multiGC_FFA_to_PFC_mean;
                allGC_toPFC_std(count) = results_NP.multiGC_FFA_to_PFC_std;
                
                colorDecoding_toPFC(count,:) = allGC_toPFC(count);
                count = count + 1;
            end
        end
    end
    
    % Extract data for plotting
    x_a = modelParams(:,1); %values for different a 
    y_c = modelParams(:,3); %values for different c
    a_Vals = reshape(x_a,[21,21])';
    c_Vals = reshape(y_c, [21,21])';
    
    decodeVals = reshape(allGC_toPFC,[21,21])';
    
    %Plot using surf plot
    figure(1)
    subplot(2,3,kk)
    surf(a_Vals,c_Vals,decodeVals); hold on; caxis([50 90]); colormap(jet)
    set(gca, 'FontSize',20);xlabel('XY_F = '); ylabel('YY_F = YY_R');
        title(['XX_F = XX_R = ' num2str(poss_d)]); zlim([45 95]);
    hold off;
    
    %Plot using imagesc
    figure(2)
    subplot(2,3,kk)
    imagesc(decodeVals); hold on; caxis([50 90]); colormap(jet)
    set(gca, 'FontSize',20, 'XTick', [1,11,21], 'XTickLabel', [poss_a(1),poss_a(11),poss_a(21)], ...
        'YTick', [1,11,21], 'YTickLabel', [poss_c(1),poss_c(11),poss_c(21)]); xlabel('XY_F = '); ylabel('YY_F = YY_R');
    title(['XX_F = XX_R = ' num2str(poss_d)]);
    hold off;
    
end