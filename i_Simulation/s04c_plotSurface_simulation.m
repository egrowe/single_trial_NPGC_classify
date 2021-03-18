%% Plot simulation classification accuracy using SURFACE plot
clear all; close all;
allPoss_d = -0.5: 0.05: 0.5; %set d outside (to allow this to be fixed at specific value)

for kk = 1:length(allPoss_d)
    
    poss_a = -0.5: 0.05: 0.5; %weight for 'a' (= x self conn X)
    poss_c = -0.5: 0.05: 0.5; %weight for 'c' (= x self conn X in Y data)
    poss_d = allPoss_d(kk) %weight for 'd' (= y self conn Y)
    
    count = 1;
    
    for gg = 1:length(poss_a)
        
        for  ff = 1:length(poss_c)
            
            for vv = 1:length(poss_d)
                
                %Model parameters
                this_a = poss_a(ss); % use this weight for 'a' (= x self conn X)
                this_c = poss_c(ff);  % use this weight for 'c' (= x self conn X in Y data)
                this_d = poss_d(vv); % use this weight for 'd' (= y self conn Y)
                
                modelParams(count,:) = [this_a, this_c, this_d]
                
                loadname = ['Results_Decode_FvR_NonPara_GC_a=' num2str(this_a) '_c=' num2str(this_c) '_d=' num2str(this_d) '.mat']
                load(loadname)
                
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
    a_Vals = [x_a(1:5)'; x_a(6:10)'; x_a(11:15)'; x_a(16:20)'; x_a(21:25)'];
    c_Vals = [y_c(1:5)'; y_c(6:10)'; y_c(11:15)'; y_c(16:20)'; y_c(21:25)'];
    
    decodeVals = [allGC_toPFC(1:5); allGC_toPFC(6:10); allGC_toPFC(11:15); allGC_toPFC(16:20); allGC_toPFC(21:25)];
    
    %Plot using surf plot
    figure(1)
    subplot(2,3,kk)
    surf(a_Vals,c_Vals,decodeVals); hold on; caxis([45 95]); colormap(jet)
    set(gca, 'FontSize',20);xlabel('cF = '); ylabel('a_F = a_R');
        title(['d_F = d_R = ' num2str(poss_d)]); zlim([45 95]);
    hold off;
    
    %Plot using imagesc
    figure(2)
    subplot(2,3,kk)
    imagesc(decodeVals); hold on; caxis([45 95]); colormap(jet)
    set(gca, 'FontSize',20, 'XTick', [1,3,5], 'XTickLabel', [poss_a(1),poss_a(3),poss_a(5)], ...
        'YTick', [1,3,5], 'YTickLabel', [poss_c(1),poss_c(3),poss_c(5)]); xlabel('c_F = '); ylabel('a_F = a_R');
    title(['d_F = d_R = ' num2str(poss_d)]);
    hold off;
    
end