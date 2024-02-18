%% This is the parameter obtain from hyperparameter tuning (From scratch)
% Define training parameters
maxEpochs = 729;                                % Maximum number of training epochs.
miniBatchSize = 11;                             % Mini batch size.
alpha0 = 0.001;                                % InitialLearnRate
tau = 10;                                      % LearnRateDropPeriod
gamma = 0.9772;                                 % LearnRateDropFactor
VBfreq = floor(size(X_train,1)/miniBatchSize);  % ValidationFrequency and VerboseFrequency
resetInputNormalization = false;                % ResetInputNormalization
gradientThreshold = 1.0906;                          % GradientThreshold
shuffle = 'every-epoch';
outputNetwork = 'best-validation-loss';
