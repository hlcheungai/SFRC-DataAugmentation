% final output data structrue:
% X:
% Row 1: a11, 2: a22, 3:a33, 4: a12, 5: a13, 6: a23 7: vf, 8: e11, 9: e22, 10: e33, 11: e12, 12: e23, 13: e13
% y: 
% Row 1: s11, 2: s22, 3:s33, 4: s12, 5: s23, 6: s13 

clear; clc;
%% Setting
name = "analysis_";
path = "C:\Users\hlche\Documents\Chalmers\Thesis\Justin_Cheung_Files\Data_Generation\FFT_Result\3D\";
nData = 250;

timeSteps = 101;        % Will use interpolation to enforce 100 even timesteps
stressUnitConv = 1e-6;  % From Pa to MPa

trainSize = 0.8;
validationSize = 0.15;
testSize = 1-trainSize-validationSize;

stressThreshold = 1e-5;
saveName = "RVE_FFT_2D_Data";
skipSave = 1;

%% Import
X = {};
Y = {};
for i=1:nData
    tmpDataName = append(path,name,int2str(i));
    try
        [vf, a11, a22, a33, a12, a13, a23] = importProp_FFT_log(append(tmpDataName,".log"), [1, Inf]);
        strainStress = importStrainStress_FFT_mac(append(tmpDataName,".mac"), [3, Inf], timeSteps);
    catch
        disp(append("Cannot find ",tmpDataName))
        continue
    end

    if size(strainStress,1)==timeSteps
        X{end+1} = [repmat([a11 a22 a33 a12 a13 a23 vf], timeSteps, 1) strainStress(:,1:6)]';
        Y{end+1} = strainStress(:,7:12)'*stressUnitConv;
    end
    progress = i/nData
end
X = X';
Y = Y';

%% Check data
maxStressHistory = [];
for i = 1:length(X)
    maxStressHistory(end+1) = max(max(abs(Y{i})));
end
X(maxStressHistory<stressThreshold)=[];
Y(maxStressHistory<stressThreshold)=[];

%% Split train, test, validation data
if trainSize+testSize >1
    disp("Error: trainSize + testSize > 1")
    return
end
nTrainData = round(length(X)*trainSize);
nTestData = round(length(X)*testSize);
nValidation = length(X)-nTrainData-nTestData;
if nValidation == 0
    disp("Warning: No validation data")
end
X_train = X(1:nTrainData);
Y_train = Y(1:nTrainData);
X_test = X(nTrainData+1:nTrainData+nTestData);
Y_test = Y(nTrainData+1:nTrainData+nTestData);
X_valid = X(nTrainData+nTestData+1:end);
Y_valid = Y(nTrainData+nTestData+1:end);

%% Save
if skipSave==0
    save(saveName, 'X_train', 'Y_train', 'X_test', 'Y_test', 'X_valid', 'Y_valid')
end