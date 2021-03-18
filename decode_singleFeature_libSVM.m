function [decodability, cv_acc] = decode_singleFeature_libSVM(data_faces, data_random, nReps, costRange, zscore)

% Classifies using SVM (libSVM) for single feature decoding
%
% Conducts 10-fold cross-validation of training data and repeats decoding
% by number of repitions specified as 'nReps'

% Conducts grid search for optimised cost parameter ('c') at each iteration
%
% Normalises data if z-score is set to 1
%
% Assumes equal number of observations in each class
%
% Inputs:
%   data_faces: 2D array (nFeatures x nObs)
%   data_random: 2D array (nFeatures x nObs)
%   nReps: number of times to repeat decoding analysis (test phase) (default: 10)
%   cost_range: vector; list of costs (actual costs, e.g. 2^-10, specify '-10') (default -20:5:5)
%   zscore: if set to 1, applies z-score transform (set to 0 for no z-score) (default zcore = 1)

% Outputs:
%   decodability: array of classification results (nFeatures x nReps)
%   cv_acostCounter: vector of cross-valisation acostCounteruracy (nFeatures x nReps)

%Randomise by clock
rand('seed', sum(100 * clock)); % seed random nuber generator by clock

% Test cost range
if ~exist('costRange','var')
    costRange = [-20:5:5]; % default cost range if not specified
end

% Set nReps if not specified
if ~exist('nReps', 'var')
    nReps = 10;
end

%Set z score to = 1 if not specified
if ~exist('zscore','var')
    zscore = 1;
end

%Gather data for random (in case we need to split if diff sizes)
data_randomOri = data_random;

for reps = 1:nReps
    
    if size(data_faces,2) ~= size(data_randomOri,2)
        data_random = data_randomOri(:,randperm(size(data_randomOri,2),size(data_faces,2)));
    end
    
    %Set trials in test/train
    trainSplit = ceil(size(data_faces,2)*0.7); %
    randSplitF = randperm(size(data_faces,2)); %select face trials to use randimly
    train70F = randSplitF(1:trainSplit); % select 70% train trials
    test30F = randSplitF(trainSplit+1:end); % select 30% test trials
    
    randSplitR = randperm(size(data_random,2)); %select random trials to use randimly
    train70R = randSplitR(1:trainSplit); % select 70% train trials
    test30R = randSplitR(trainSplit+1:end); % select 30% test trials
    
    %Get test and train sets
    trainF = data_faces(:,train70F); % features x trials
    trainR = data_random(:,train70R);% features x trials
    testF = data_faces(:,test30F);% features x trials
    testR = data_random(:,test30R);% features x trials
    
    % %Combine the parameters to z-score TRAINING SET first across face and non-face trials
    allTrain = [trainF,trainR]; %features x trials (F/R train)
    allTest = [testF, testR]; % features x trials (F/R test)
%     allData_orig = [trainF, testF, trainR, testR]; % features x trials (all)
%     allData_Train_Test = [allTrain, allTest];
%     
    %Add labels for each row of features to be decoded (single feature decoding)
    i = 1;
    for iCon = 1:size(data_faces,1)
        classLabel_train(i,:) = [repmat(iCon, 1, size(trainF,2)), repmat(iCon*10, 1, size(trainR,2))] ;
        i = i + 1;
    end
    
    i = 1;
    for iCon = 1:size(data_faces,1)
        classLabel_test(i,:) = [repmat(iCon, 1, size(testF,2)), repmat(iCon*10, 1, size(testR,2))] ;
        i = i + 1;
    end
    
    if zscore == 1
        %% Zscore from here
        %Zscore across each row for TRAINING SET (all combined)
%         for ii = 1:size(allTrain,1)
%             Xtrain = allTrain(ii,:); %allTrain = features x nTrials
%             zXtrain(ii,:) = (Xtrain-mean(Xtrain))./std(Xtrain);
%         end
%         
                    % Normalise training data
            means = mean(allTrain, 2);
            stds = std(allTrain, [], 2);
            means_mat = repmat(means, 1,[size(allTrain, 2)]);
            stds_mat = repmat(stds, 1,[size(allTrain, 2)]);
            allTrain = (allTrain - means_mat) ./ stds_mat;
        
%         %Zscore across each row for TESTING SET (all combined)
%         for jj = 1:size(allTest,1)
%             Xtrain = allTrain(jj,:);
%             Ytest = allTest(jj,:);
%             zYTest(jj,:) = (Ytest-mean(Xtrain))/std(Xtrain);
%         end
        
                    % Normalise validation data (using same parameters as for
            % training data)
            means_mat = repmat(means, 1, [size(allTest, 2)]);
            stds_mat = repmat(stds, 1, [size(allTest, 2)]);
            allTest = (allTest - means_mat) ./ stds_mat;
        
        %Re-divide the data by face train/test and random train/test
        zfTrain = allTrain(:,1:size(trainF,2)); %zscored faces training data
        zrTrain = allTrain(:,size(trainF,2)+1:end); %zscored random training data
        zfTest = allTest(:,1:size(testF,2)); %zscore faces test data
        zrTest = allTest(:,size(testF,2)+1:end); %zscored random test data
        
        clear allTrain allTest
        allTrain = [zfTrain, zrTrain];
        allTest = [zfTest, zrTest];
    end
    
    %% Optimise the cost parameter through grid search
    for freqBand = 1:size(classLabel_train,1)
        thisFreqCV(freqBand,:) = 0;
        for costCounter = 1:length(costRange) %vary cost expoentially
            log2c = costRange(costCounter);
            cmdCV = ['-s 0 -t 0 -c ' num2str(2^(log2c)) ' -v 10 -q']; %
            
            labels = classLabel_train(freqBand,:)';
            data = allTrain(freqBand,:)';
            
            %CROSS-VALIDATE THE CLASSIFIER (10-fold)
            [cv] = svmtrain(labels, data , cmdCV);
            cv_acc(freqBand,costCounter,reps) = cv;
            
            if (cv >= thisFreqCV(freqBand,:))
                thisFreqCV(freqBand,:) = cv; bestcv(freqBand,:) = log2c;% bestg = 2^log2g;
            end
            fprintf('%g %g (best c=%g,)\n', log2c, thisFreqCV(freqBand,:), bestcv(freqBand,:));
        end
        allBestCostParams(freqBand,reps) = bestcv(freqBand,:);
    end
    
    cParamBest(reps) = mode(allBestCostParams(:,reps)); %selec the most frequent, BEST CV across channels (use for decoding)
    cIdx = costRange == cParamBest(reps); % find index for best CV
    useCIdx = find(cIdx == 1);
    
    %% TRAIN USING BEST COST PARAM
    cmdModel = ['-s 0 -t 0 -c ' num2str(2^costRange(useCIdx)) ' -q']; %
   
    for freqBand = 1:size(classLabel_train,1)
        labels = classLabel_train(freqBand,:)';
        data = allTrain(freqBand,:)';
        
        %TRAIN THE CLASSIFIER
        model(freqBand) = svmtrain(labels, data , cmdModel);
    end
    

    %% DECODE USING BEST COST PARAM
    %TEST THE CLASSIFIER
    for freqBand = 1:size(allTest,1)
        [predicted_label, accuracy, decision_values]  = svmpredict(classLabel_test(freqBand,:)',allTest(freqBand,:)',model(:,freqBand));
        all_accuracy(freqBand,:) = accuracy(1,1)
         clear predicted_label accuracy decision_values
    end
    
    decodability(:,reps) = all_accuracy;
    
    clearvars -except allBestCostParams data_randomOri visualise decodability crossValAcc cv data_faces ...
        data_random cParamBest costRange accuracy cv_acc length70Split classLabel reps nReps data_faces ...
        data_random zscore tStatResults useCIdx
end
end