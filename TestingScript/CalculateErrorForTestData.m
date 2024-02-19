clear; clc; 
%% Load network
load("FromScratch_RR_15_best.mat")

rotate2PrincipalAxis = 0;
timeSteps = "original";

%% Load data
ALL = load("RVE_all_data.mat","X_test","Y_test");

% U = load("GeneralTest_Uniaxial.mat","X_test","Y_test");
% S = load("GeneralTest_Shear.mat","X_test","Y_test");
% PS = load("GeneralTest_PlainStrain.mat","X_test","Y_test");
% NS = load("GeneralTest_BiNormalShear.mat","X_test","Y_test");
% NN = load("GeneralTest_BiNormalNormal.mat","X_test","Y_test");
% ALL.X_test = [U.X_test; S.X_test; PS.X_test; NS.X_test; NN.X_test];
% ALL.Y_test = [U.Y_test; S.Y_test; PS.Y_test; NS.Y_test; NN.Y_test];

%% Prediction and calculate error
X_test = ALL.X_test;
Y_test = ALL.Y_test;
vonMisesMeRE_History = [];
vonMisesMaRE_History = [];
for j = 1:length(X_test)
    input = X_test{j};
    target = Y_test{j};
    switch timeSteps
        case "original"
        otherwise
            tmpTimeSteps = str2double(timeSteps);
            input = interp1(linspace(0,1,size(input,2)),input',linspace(0,1,tmpTimeSteps))';
    end

    [inputForNN, ~, r] = preformPrincipalRotation(input,target,rotate2PrincipalAxis);
    predictionFromNN = predict(net,inputForNN);
    [~, prediction] = rotateBackToOrginalCoordinate(inputForNN,predictionFromNN,r);
    prediction = interp1(linspace(0,1,size(prediction,2)), prediction', linspace(0,1,size(target,2)))';

    vonMisesTarget = calculateVonMisesStress(target);
    vonMisesPrediction = calculateVonMisesStress(prediction);
    vonMisesError = vonMisesPrediction - vonMisesTarget;
    vonMisesMeRE = sqrt(sum(vonMisesError.^2,2)/length(vonMisesError))/range(vonMisesTarget);
    vonMisesMaRE = max(abs(vonMisesError))/range(vonMisesTarget);

    vonMisesMeRE_History(end+1) = vonMisesMeRE;
    vonMisesMaRE_History(end+1) = vonMisesMaRE;
end
meanVonMisesMeRE = mean(vonMisesMeRE_History);
meanVonMisesMaRE = mean(vonMisesMaRE_History);

%% function
function vonMisesStress = calculateVonMisesStress(stress)
s11 = stress(1,:);
s22 = stress(2,:);
s33 = stress(3,:);
s12 = stress(4,:);
s23 = stress(5,:);
s13 = stress(6,:);
vonMisesStress = sqrt(((s11-s22).^2 + (s22-s33).^2 + (s33-s11).^2 + 6.*(s12.^2+s23.^2+s13.^2))/2);
end

function [xRotated, yRotated, r] = preformPrincipalRotation(x,y,confirmRotation)
if confirmRotation > 0
    timeSteps = size(x,2);
    vf = x(7,1);

    a11 = x(1,1); a22 = x(2,1);  a33 = x(3,1);
    a12 = x(4,1); a13 = x(5,1);  a23 = x(6,1);
    tmpA = [a11 a12 a13; a12 a22 a23; a13 a23 a33];
    [tmpR,~] = eig(tmpA);
    r = tmpR';

    ARotated=r*tmpA*r';

    e11 = x(8,:);  e22 = x(9,:);  e33 = x(10,:);
    e12 = x(11,:)/2; e23 = x(12,:)/2; e13 = x(13,:)/2;
    tmpStrain = [e11; e12; e13; e12; e22; e23; e13; e23; e33];
    tmpStrain = reshape(tmpStrain,[3,3,timeSteps]);

    s11 = y(1,:); s22 = y(2,:); s33 = y(3,:);
    s12 = y(4,:); s23 = y(5,:); s13 = y(6,:);
    tmpStress = [s11; s12; s13; s12; s22; s23; s13; s23; s33];
    tmpStress = reshape(tmpStress,[3,3,timeSteps]);

    r = repmat(r,1,1,timeSteps);

    strainRotated = pagemtimes(pagemtimes(r,tmpStrain),'none',r,'transpose');
    strainRotated = reshape(strainRotated,[9,size(strainRotated,3)])';

    stressRotated = pagemtimes(pagemtimes(r,tmpStress),'none',r,'transpose');
    stressRotated = reshape(stressRotated,[9,size(stressRotated,3)])';

    prop = [ARotated(1,1) ARotated(2,2) ARotated(3,3) ARotated(1,2) ARotated(1,3) ARotated(2,3) vf];

    strainStress = [strainRotated(:,1), strainRotated(:,5), strainRotated(:,9), strainRotated(:,2)*2, strainRotated(:,6)*2, strainRotated(:,3)*2,...
    stressRotated(:,1), stressRotated(:,5), stressRotated(:,9), stressRotated(:,2), stressRotated(:,6), stressRotated(:,3)];
    
    xRotated = [repmat(prop, timeSteps, 1) strainStress(:,1:6)]';
    yRotated = strainStress(:,7:12)';
    r = r(:,:,1);
else
    xRotated = x;
    yRotated = y;
    r = [];
end
end

function [xRotated, yRotated] = rotateBackToOrginalCoordinate(x,y,r)
if isempty(r)
    xRotated = x;
    yRotated = y;
else
    r = r';

    timeSteps = size(x,2);
    vf = x(7,1);

    a11 = x(1,1); a22 = x(2,1);  a33 = x(3,1);
    a12 = x(4,1); a13 = x(5,1);  a23 = x(6,1);
    tmpA = [a11 a12 a13; a12 a22 a23; a13 a23 a33];
    
    ARotated=r*tmpA*r';

    e11 = x(8,:);  e22 = x(9,:);  e33 = x(10,:);
    e12 = x(11,:)/2; e23 = x(12,:)/2; e13 = x(13,:)/2;
    tmpStrain = [e11; e12; e13; e12; e22; e23; e13; e23; e33];
    tmpStrain = reshape(tmpStrain,[3,3,timeSteps]);

    s11 = y(1,:); s22 = y(2,:); s33 = y(3,:);
    s12 = y(4,:); s23 = y(5,:); s13 = y(6,:);
    tmpStress = [s11; s12; s13; s12; s22; s23; s13; s23; s33];
    tmpStress = reshape(tmpStress,[3,3,timeSteps]);

    r = repmat(r,1,1,timeSteps);

    strainRotated = pagemtimes(pagemtimes(r,tmpStrain),'none',r,'transpose');
    strainRotated = reshape(strainRotated,[9,size(strainRotated,3)])';

    stressRotated = pagemtimes(pagemtimes(r,tmpStress),'none',r,'transpose');
    stressRotated = reshape(stressRotated,[9,size(stressRotated,3)])';

    prop = [ARotated(1,1) ARotated(2,2) ARotated(3,3) ARotated(1,2) ARotated(1,3) ARotated(2,3) vf];

    strainStress = [strainRotated(:,1), strainRotated(:,5), strainRotated(:,9), strainRotated(:,2)*2, strainRotated(:,6)*2, strainRotated(:,3)*2,...
    stressRotated(:,1), stressRotated(:,5), stressRotated(:,9), stressRotated(:,2), stressRotated(:,6), stressRotated(:,3)];
    
    xRotated = [repmat(prop, timeSteps, 1) strainStress(:,1:6)]';
    yRotated = strainStress(:,7:12)';
end
end
