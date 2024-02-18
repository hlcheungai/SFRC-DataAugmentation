clear; clc;
%% Define optimization parameter and range
maxEpochs = optimizableVariable('maxEpochs',[300,1000],'Type','integer');               % Maximum number of training epochs
miniBatchSize = optimizableVariable('miniBatchSize',[5,100],'Type','integer');       % Mini batch size
alpha0 = optimizableVariable('alpha0',[0.0001,0.001],'Type','real');                        % InitialLearnRate
tau = optimizableVariable('tau',[10,200],'Type','integer');                           % LearnRateDropPeriod
gamma = optimizableVariable('gamma',[0.8,0.99],'Type','real');                          % LearnRateDropFactor

%% Settings
dataset1 = "RVE_all_data.mat"
dataset2 = "RVE_all_data_RR_9.mat"
% display = 'training-progress';
display = 'none';

nTrial = 30

folderPath = fullfile(pwd,"Experiment_RR_9_second")
%folderPath = pwd;
baseFileName = 'trial';

%% Optimization
FinalValidationLoss = @(params) trainNetworkAndReturnValidationLoss(params,dataset1,dataset2,folderPath,baseFileName,display);
optimizationResults = bayesopt( ...
    FinalValidationLoss, ...
    [maxEpochs,miniBatchSize,alpha0,tau,gamma], ...
    "MaxObjectiveEvaluations",nTrial);

%% Rename best network
[~,bestIteration] = min(optimizationResults.ObjectiveMinimumTrace);
oldName = fullfile(folderPath, append("trial_",int2str(bestIteration),".mat"));
newName = fullfile(folderPath, append("trial_",int2str(bestIteration),"_best.mat"));
movefile(oldName, newName);

%% Save settings
tempPath = fullfile(folderPath, "experiment_setting.mat");
save(tempPath,"optimizationResults")
