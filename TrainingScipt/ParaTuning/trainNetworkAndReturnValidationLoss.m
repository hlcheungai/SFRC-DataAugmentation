function FinalValidationLoss = trainNetworkAndReturnValidationLoss(params,dataset1,dataset2,folderPath,baseFileName,display)
%% Import hyper-parameters to try
maxEpochs = params.maxEpochs;                   % Maximum number of training epochs.
miniBatchSize = params.miniBatchSize;           % Mini batch size.
alpha0 = params.alpha0;                         % InitialLearnRate
tau = params.tau;                               % LearnRateDropPeriod
gamma = params.gamma;                           % LearnRateDropFactor
gradientThreshold = 1;                          % GradientThreshold

%% Import data
original = load(dataset1);
rotated = load(dataset2);
X_train = [original.X_train; rotated.X_train];
Y_train = [original.Y_train; rotated.Y_train];
X_valid = [original.X_valid; rotated.X_valid];
Y_valid = [original.Y_valid; rotated.Y_valid];

%% Define training options
VBfreq = floor(size(X_train,1)/miniBatchSize);
options = trainingOptions('adam', ...
    'ExecutionEnvironment','auto', ...
    'MaxEpochs',maxEpochs, ...
    'Shuffle','every-epoch', ...
    'ResetInputNormalization', false,...
    'MiniBatchSize',miniBatchSize, ...
    'GradientThreshold',gradientThreshold, ...
    'InitialLearnRate',alpha0, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',tau, ...
    'LearnRateDropFactor',gamma, ...
    'Verbose',true, ...
    'VerboseFrequency',VBfreq, ...
    'Plots',display, ...
    'ValidationFrequency',VBfreq, ...
    'ValidationData',{X_valid,Y_valid}, ...
    'OutputNetwork','best-validation-loss');

layers = [
    sequenceInputLayer(13,"Name","sequence","Normalization","zscore")
    gruLayer(500,'Name','gru_1')
    gruLayer(500,'Name','gru_2')
    gruLayer(500,'Name','gru_3')
    dropoutLayer(0.5,'Name','dropout')
    fullyConnectedLayer(6,'Name','fc')
    regressionLayer("Name","regressionoutput")];
%% Train network and get final validation loss
[net,info] = trainNetwork(X_train,Y_train,layers,options);
FinalValidationLoss = info.FinalValidationLoss;

%% Save file
index = 1;
while true
    fileName = sprintf('%s_%d.mat', baseFileName, index);
    if exist(fullfile(folderPath, fileName), 'file')
        index = index + 1;
    else
        break;
    end
end
save(fullfile(folderPath, fileName), "net","info","params","dataset1","dataset2");
