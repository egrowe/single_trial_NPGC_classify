function [decodability, cv_acc, w] = decode_multiFeature_libSVM(data_faces, data_random, nReps, costRange, zVal)

% Classifies using SVM (libSVM) for multifeature datasets
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
%   nReps: nunmber of times to repeat decoding analysis (test phase) (default: 10)
%   cost_range: vector; list of costs (actual costs, e.g. 2^-10, specify '-10') (default -20:5:5)
%   zscore: if set to 1, applies z-score transform (set to 0 for no z-score) (default zcore = 1)

% Outputs:
%   decodability: vector of classification results (cost_range x nReps)
%   cv_acc: vector of cross-valisation accuracy (cost_range x nReps)


% Test cost range
if ~exist('costRange','var')
    costRange = [-20:5:5]; % default cost range if not specified
end

% Set nReps if not specified
if ~exist('nReps', 'var')
    nReps = 10;
end

%Set z score to = 1 if not specified
if ~exist('zVal','var')
    zVal = 1;
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
    
    % %Add labels for each row of features to be decoded (single featuredecoding)
    classLabel_train = [repmat(1, 1, size(trainF,2)), repmat(10, 1, size(trainR,2))]';
    classLabel_test = [repmat(1, 1, size(testF,2)), repmat(10, 1, size(testR,2))]' ;
    
    if zVal == 1
        %% Zscore normalise from here
        % Normalise training data
        means = mean(allTrain, 2);
        stds = std(allTrain, [], 2);
        means_mat = repmat(means, 1,[size(allTrain, 2)]);
        stds_mat = repmat(stds, 1,[size(allTrain, 2)]);
        allTrain = (allTrain - means_mat) ./ stds_mat;
        
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
        
        allZTrain = [zfTrain, zrTrain];
        allZTest = [zfTest, zrTest];
        clear allTrain allTest
        
        % Use Z-score data if specified at beginnig of script
        allTrain = allZTrain;
        allTest = allZTest;
    end
    
    %% Optimise the cost parameter through grid search
    % bestcv = [];
    thisChanCV = 0;
    for costCounter = 1:length(costRange)
        log2c = costRange(costCounter);
        
        cmdCV = ['-s 0 -t 0 -c ' num2str(2^log2c) ' -v 10 -q']; %
        labels = classLabel_train;
        data = allTrain';
        
        %CROSS-VALIDATE THE CLASSIFIER (10-fold)
        cv = svmtrain(labels, data , cmdCV);
        cv_acc(reps,costCounter) = cv;
        
        if (cv >= thisChanCV)
            thisChanCV = cv; bestcv = log2c;% bestg = 2^log2g;
        end
        % fprintf('%g %g (best c=%g,)\n', log2c, thisChanCV(reps,:), bestcv(reps,:));
    end
    allBestCostParams(reps) = bestcv;
    
    cParamBest(reps) = mode(allBestCostParams(:,reps)); %selec the most frequent, BEST CV across channels (use for decoding)
    cIdx = costRange == cParamBest(reps); % find index for best CV
    useCIdx = find(cIdx == 1);

    %% TRAIN THE CLASSIFIER
    cmdModel = ['-s 0 -t 0 -c ' num2str(2^costRange(useCIdx)) ' -q']; %
    %TRAIN THE CLASSIFIER
    model(reps) = svmtrain(labels, data , cmdModel);
    w(reps,:) = model(reps).SVs' * model(reps).sv_coef;
    
  
    %% TEST THE CLASSIFIER
        [predicted_label, accuracy, decision_values]   = svmpredict(classLabel_test,allTest',model(reps));
        decodability(reps) = accuracy(1,:);
    
    clearvars -except w data_randomOri decodability allBestCostParams accuracy crossValAcc allZTrain allZTest useCParam ...
        data_faces allTest allTrain classLabel_train classLabel_test allTrain useCIdx ...
        data_random costRange cv_acc length70Split classLabel rr nReps data_faces data_random zVal
end
end