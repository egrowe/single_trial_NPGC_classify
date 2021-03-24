function [decodability, cv_acc] = decode_singleFeature_libSVM(data_faces, data_random, nReps, costRange, zVal)

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
%   data_faces: 2D array (nObs x nFeatures)
%   data_random: 2D array (nObs x nFeatures)
%   nReps: number of times to repeat decoding analysis (test phase) (default: 10)
%   cost_range: vector; list of costs (actual costs, e.g. 2^-10, specify '-10') (default -20:5:5)
%   zscore: if set to 1, applies z-score transform (set to 0 for no z-score) (default zcore = 1)

% Outputs:
%   decodability: array of classification results (nFeatures x nReps)
%   cv_acostCounter: vector of cross-valisation acostCounteruracy (nFeatures x nReps)

%Check if data input is correct
if size(data_faces,2) > size(data_faces,1) 
    error('Error! Input data incorrect format')
end

%Training split proportion (here 70%)
trainProp = 0.7; % percentage to split as training set

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
    zVal = 1;
end

%Gather data for random (in case we need to split if diff sizes)
data_random_all = data_random;

for reps = 1:nReps
    
    %Randomise the random number generator each loop
    rng('shuffle')
    
    if size(data_faces,1) ~= size(data_random_all,1)
        data_random = data_random_all(randperm(size(data_random_all,1),size(data_faces,1)));
    end
    
    %Set trials in test/train
    trainSplit = ceil(size(data_faces,1)*trainProp); %number of trials for 70%
    
    %Partition the FACE data randomly (randomise trial order)
    randTrials_F = randperm(size(data_faces,1));  %select trial order randomly 
    trainF = data_faces(randTrials_F(1:trainSplit),:); % select 70% train trials: Obs x Features
    testF = data_faces(randTrials_F(trainSplit+1:end),:); % select 30% test trials: Obs x Features
    
   %Partition the RANDOM data randomly (randmomise trial order)
    randTrials_R = randperm(size(data_random,1));  %select trial order randomly
    trainR = data_random(randTrials_R(1:trainSplit),:) % select 70% train trials: Obs x Features
    testR = data_random(randTrials_R(trainSplit+1:end),:); % select 30% test trials: Obs x Features

    % Gather all TRAIN and TEST data
    allTrain = [trainF;trainR]; %features x trials (F/NF train)
    allTest = [testF; testR]; % features x trials (F/NF test)

    % Gather all LABELS for TRAIN and TEST sets
   classLabel_train = [repmat(0, size(trainF,1),1); repmat(1, size(trainR,1),1)];
   classLabel_test = [repmat(0, size(testF,1),1); repmat(1, size(testR,1),1)] ;
    
    if zVal == 1
            % Normalise training data
            means = mean(allTrain, 1);
            stds = std(allTrain, [], 1);
            means_mat = repmat(means, [size(allTrain, 1), 1]);
            stds_mat = repmat(stds, [size(allTrain, 1), 1]);
            allTrain = (allTrain - means_mat) ./ stds_mat;
            
            % Normalise validation data (using same parameters as for
            % training data)
            means_mat = repmat(means, [size(allTest, 1), 1]);
            stds_mat = repmat(stds, [size(allTest, 1), 1]);
            allTest = (allTest - means_mat) ./ stds_mat;
    end
    
    %% Optimise the cost parameter through grid search
    for freqBand = 1:size(allTrain,2)
        thisFreqCV(freqBand,:) = 0;
        for costCounter = 1:length(costRange) %vary cost expoentially
            cVal = costRange(costCounter);
            log2c = log(2)^cVal;

            cmdCV = ['-s 0 -t 0 -c ' num2str(log2c) ' -v 10 -q']; %
            
            labels = classLabel_train;
            data = allTrain(:,freqBand);
            
            %CROSS-VALIDATE THE CLASSIFIER (10-fold)
            [cv] = svmtrain(labels, data , cmdCV);
            cv_acc(freqBand,costCounter,reps) = cv;
            
            if (cv >= thisFreqCV(freqBand,:))
                thisFreqCV(freqBand,:) = cv; bestcv(freqBand,:) = cVal;% bestg = 2^log2g;
            end
            fprintf('%g %g (best c=%g,)\n', cVal, thisFreqCV(freqBand,:), bestcv(freqBand,:));
        end
        allBestCostParams(freqBand,reps) = bestcv(freqBand,:);
    end
    
    cParamBest(reps) = mode(allBestCostParams(:,reps)); %selec the most frequent, BEST CV across channels (use for decoding)

    %% TRAIN USING BEST COST PARAM
   
    for freqBand = 1:size(allTrain,2)
        %Select best C parameter for this frequency
        useLog2c = log(2)^allBestCostParams(freqBand); %choose best C parameter
        cmdModel = ['-s 0 -t 0 -c ' num2str(useLog2c) ' -q']; %

        %Train the model using this parameter
        labels = classLabel_train;
        data = allTrain(:,freqBand);
        
        %TRAIN THE CLASSIFIER
        trained(freqBand) = svmtrain(labels, data , cmdModel);
    end
    

    %% DECODE USING BEST COST PARAM
    %TEST THE CLASSIFIER
    for freqBand = 1:size(allTest,2)
        [predicted_label, accuracy, decision_values]  = svmpredict(classLabel_test,allTest(:,freqBand),trained(1,freqBand));
        all_accuracy(freqBand,:) = accuracy(1,1)
         clear predicted_label accuracy decision_values
    end
    
    decodability(:,reps) = all_accuracy;
    
     clearvars -except w data_randomOri decodability allBestCostParams crossValAcc useCParam ...
        data_faces data_random_all trainProp data_random costRange zVal model nReps cv_acc trained
end
end