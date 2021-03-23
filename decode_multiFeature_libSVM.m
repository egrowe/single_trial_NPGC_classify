function [decodability, cv_acc, w] = decode_multiFeature_libSVM_new(data_faces, data_random, nReps, costRange, zVal)

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
%   data_faces: 2D array (nObs n nFeatures)
%   data_random: 2D array (nObs n nFeatures)
%   nReps: nunmber of times to repeat decoding analysis (test phase) (default: 10)
%   cost_range: vector; list of costs (actual costs, e.g. 2^-10, specify '-10') (default -20:5:5)
%   zscore: if set to 1, applies z-score transform (set to 0 for no z-score) (default zcore = 1)

% Outputs:
%   decodability: vector of classification results (cost_range x nReps)
%   cv_acc: vector of cross-valisation accuracy (cost_range x nReps)

%Training split proportion (here 70%)
trainProp = 0.70; % percentage to split as training set

% %Check if data input is correct
% if size(data_faces,2) > size(data_faces,1) 
%     error('Error! Input data incorrect format')
% end

% Test cost range
if ~exist('costRange','var')
    costRange = [-15:5:15]; % default cost range if not specified
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
data_random_all = data_random;

for reps = 1:nReps
    
    %Reset the random number generator by the internal clock each loop
    t = clock;
    rng(t(3) * t(4) * t(5),'twister')
    
        %Shuffle datasets
    data_faces = Shuffle(data_faces'); data_faces = data_faces';
    data_random_all = Shuffle(data_random_all'); data_random_all = data_random_all';
    data_random = Shuffle(data_random'); data_random = data_random';
    
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
    % bestcv = [];
    thisChanCV = 0;
    for costCounter = 1:length(costRange)
        cVal = costRange(costCounter);
        log2c = log(2)^cVal;
        
        %Setup the model parameters and define data + labels
        cmdCV = ['-s 0 -t 0 -c ' num2str(log2c) ' -v 10 -q']; %
        labels = classLabel_train;
        data = allTrain;
        
        % 10-fold CROSS VALIDATION TO TUNE COST PARAMETER
        cv = svmtrain(labels, data , cmdCV); %CV accuracy
        cv_acc(reps,costCounter) = cv; %save the CV accuracy
        
        % Determine the best cost parameter (i.e. highest CV accuracy)
        if (cv >= thisChanCV)
            thisChanCV = cv; bestcv = cVal;
        end
    end
    allBestCostParams(reps) = bestcv; %Save each of the best cost parameters (per repetition of outer loop)
    
    cParamBest(reps) = mode(allBestCostParams(:,reps)); %selec the most frequent, BEST CV across channels (use for decoding)
    cIdx = costRange == cParamBest(reps); % find index for best CV
    useCIdx = find(cIdx == 1);
    useLog2c = log(2)^costRange(useCIdx); %choose best C parameter
    
    %% TRAIN THE CLASSIFIER
    cmdModel = ['-s 0 -t 0 -c ' num2str(useLog2c) ' -q']; %
    %TRAIN THE CLASSIFIER
    model(reps) = svmtrain(labels, data , cmdModel);
    thisModel = model(reps);
    w(reps,:) = model(reps).SVs' * model(reps).sv_coef;
    
  
    %% TEST THE CLASSIFIER
        [predicted_label, accuracy, decision_values] = svmpredict(classLabel_test,allTest,thisModel);
        decodability(reps) = accuracy(1,:);
    
    clearvars -except w data_randomOri decodability allBestCostParams crossValAcc useCParam ...
        data_faces data_random_all trainProp data_random costRange zVal model nReps cv_acc
end
end