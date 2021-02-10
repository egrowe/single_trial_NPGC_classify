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
    randSplitF = randperm(size(data_faces,2)); %select trials to use randimly
    train70F = randSplitF(1:trainSplit); % select 70% train trials
    test30F = randSplitF(trainSplit+1:end); % select 30% test trials
    randSplitR = randperm(size(data_random,2)); %select trials to use randimly
    train70R = randSplitR(1:trainSplit); % select 70% train trials
    test30R = randSplitR(trainSplit+1:end); % select 30% test trials
    
    %Get test and train sets
    trainF = data_faces(:,train70F); % features x trials
    trainR = data_random(:,train70R);% features x trials
    testF = data_faces(:,test30F);% features x trials
    testR = data_random(:,test30R);% features x trials
    
    % %Combine the parameters to z-score TRAINING SET first across face and non-face trials
    allTrain = [trainF,trainR]; %features x trials (F/NF train)
    allTest = [testF, testR]; % features x trials (F/NF test)
    allData_orig = [trainF, testF, trainR, testR]; % features x trials (all)
    allData_Train_Test = [allTrain, allTest];
    
    %Add labels for each row of features to be decoded (single featuredecoding)
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
    
    %% Zscore from here
    %Zscore across each row for TRAINING SET (all combined)
    for ii = 1:size(allTrain,1)
        Xtrain = allTrain(ii,:); %allTrain = features x nTrials
        zXtrain(ii,:) = (Xtrain-mean(Xtrain))./std(Xtrain);
    end
    
    %Zscore across each row for TESTING SET (all combined)
    for jj = 1:size(allTest,1)
        Xtrain = allTrain(jj,:);
        Ytest = allTest(jj,:);
        zYTest(jj,:) = (Ytest-mean(Xtrain))/std(Xtrain);
    end
    
    %Re-divide the data by face train/test and random train/test
    zfTrain = zXtrain(:,1:size(trainF,2)); %zscored faces training data
    zrTrain = zXtrain(:,size(trainF,2)+1:end); %zscored random training data
    zfTest = zYTest(:,1:size(testF,2)); %zscore faces test data
    zrTest = zYTest(:,size(testF,2)+1:end); %zscored random test data
    
    allZTrain = [zfTrain, zrTrain];
    allZTest = [zfTest, zrTest];
    
    % Use Z-score data if specified at beginnig of script
    if zscore == 1
        allTrain = allZTrain;
        allTest = allZTest;
    end
    
    %% Optimise the cost parameter through grid search
    for chan = 1:size(classLabel_train,1)
        thisChanCV(chan,:) = 0;
        for costCounter = 1:length(costRange) %vary cost expoentially
            log2c = costRange(costCounter);
            cmdCV = ['-s 0 -t 0 -c ' num2str(2^(log2c)) ' -v 10 -q']; %
            cmdModel = ['-s 0 -t 0 -c ' num2str(2^log2c) ' -q']; %
            
            labels = classLabel_train(chan,:)';
            data = allTrain(chan,:)';
            
            %CROSS-VALIDATE THE CLASSIFIER (10-fold)
            cv = svmtrain(labels, data , cmdCV);
            cv_acc(chan,costCounter,reps) = cv;
            
            %TRAIN THE CLASSIFIER
            model(chan,costCounter,reps) = svmtrain(labels, data , cmdModel);
            
            if (cv >= thisChanCV(chan,:))
                thisChanCV(chan,:) = cv; bestcv(chan,:) = log2c;% bestg = 2^log2g;
            end
            fprintf('%g %g (best c=%g,)\n', log2c, thisChanCV(chan,:), bestcv(chan,:));
        end
        allBestCostParams(chan,reps) = bestcv(chan,:);
    end
    
    cParamBest(reps) = mode(allBestCostParams(:,reps)); %selec the most frequent, BEST CV across channels (use for decoding)
    cIdx = costRange == cParamBest(reps); % find index for best CV
    useCIdx = find(cIdx == 1);

    %% DECODE USING BEST COST PARAM
    
    %TEST THE CLASSIFIER
    for chan = 1:size(allTest,1)
        [predicted_label(chan,:,reps), accuracy(chan,:,reps), decision_values(chan,:,reps)]  = svmpredict(classLabel_test(chan,:)',allTest(chan,:)',model(chan,useCIdx,reps));
    end
    
    decodability(:,reps) = accuracy(:,1,reps);
    
    clearvars -except allBestCostParams data_randomOri visualise decodability crossValAcc cv data_faces ...
        data_random cParamBest costRange accuracy cv_acc length70Split classLabel reps nReps data_faces ...
        data_random zscore tStatResults useCIdx
end
end