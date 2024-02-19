%%
% clear; clc; 
% load("RVE_all_data.mat")

% load("GeneralTest_Uniaxial.mat")
% load("GeneralTest_Shear.mat")
% load("GeneralTest_PlainStrain.mat")
% load("GeneralTest_BiNormalShear.mat")
load("GeneralTest_BiNormalNormal.mat")

% load("RVE_all_data.mat")

networkToTest = ["FromScratch_RR_15_best.mat"];
nTimeSteps = ["orginal"];
nameSetting = ["R15"];

rotate2PrincipalAxis = [0, 0, 0, 0, 0, 0];
nInputs = [13, 13, 13, 13, 13, 13];
lineSetting = ["-","-","-","-","-","-"];
colourSetting = ["#0072BD","#77AC30","r","c","g","y"];
fontSizeSetting = 10;
fontNameSetting = 'Times New Roman';

targetLineType = "none";
targetMarker = "x";
targetColour = "k";

testCase = "test"; % "train", "test", "valid";
nData = 5;
finalTime = 1;

%% Load or generate testing input
input = cell(length(networkToTest),1);
target = cell(length(networkToTest),1);
t = cell(length(networkToTest),1);
vmTarget = cell(length(networkToTest),1);
vmPrediction = cell(length(networkToTest),1);

switch testCase
    case {"train", "test", "valid"}
        if testCase == "test"
            x = X_test;
            y = Y_test;
        elseif testCase == "train"
            x = X_train;
            y = Y_train;
        elseif testCase == "valid"
            x = X_valid;
            y = Y_valid;
        end
        try
            orginalInput = x{nData};
            orginalTarget = y{nData};
            endStep = min(round(finalTime/1*size(orginalInput,2),0)+1,size(orginalInput,2));
            orginalInput = orginalInput(:,1:endStep);
            orginalTarget = orginalTarget(:,1:endStep);
        catch
            disp("Error: Invaid test number")
        end
    otherwise
        disp("Error: Invaid test case")
end
%% Preform prediction
for i = 1:length(networkToTest)
    switch nTimeSteps(i)
        case "orginal"
            tmpInput = orginalInput;
            tmpTarget = orginalTarget;
            tmpT = linspace(0,1,size(tmpInput,2));
        otherwise
            tmpInput = orginalInput;
            tmpTarget = orginalTarget;
            tmpT = linspace(0,1,size(tmpInput,2));
            tmpInput = interp1(linspace(0,1,size(tmpInput,2)), tmpInput', linspace(0,1,str2double(nTimeSteps(i))))';
    end

    [tmpInput, ~, r] = preformPrincipalRotation(tmpInput, tmpTarget,rotate2PrincipalAxis(i), nInputs(i), 1);
    input{i} = tmpInput;
    target{i} = tmpTarget;
    vmTarget{i} = calculateVonMisesStress(tmpTarget);
    t{i} = tmpT;
end     
%% Plot
plottingInput(linspace(0,1,size(orginalInput,2)), orginalInput, orginalTarget, targetColour, targetLineType, targetMarker,fontSizeSetting,fontNameSetting)
vmMeRE_History = zeros(length(networkToTest),1);
vmMaRE_History = zeros(length(networkToTest),1);
for i = 1:length(networkToTest)
    load(networkToTest(i),"net")
    tic
    prediction = predict(net,input{i});
    toc
    if rotate2PrincipalAxis(i) == 1
        prediction = rotateBackToOrginalCoordinateFromPrincipal(prediction,r);
    end
    prediction = interp1(linspace(0,1,size(prediction,2)), prediction', t{i})';
    vmPrediction{i} = calculateVonMisesStress(prediction);
    vmError = vmPrediction{i} - vmTarget{i};
    vmMeRE = sqrt(sum(vmError.^2,2)/length(vmError))/range(vmTarget{i});
    vmMaRE = max(abs(vmError))/range(vmTarget{i});
    vmPercentageError = abs(vmError./vmTarget{i})*100;

    vmMeRE_History(i) = vmMeRE;
    vmMaRE_History(i) = vmMaRE;

    subplot(2,4,1)
    plot(t{i},prediction(1,:),"LineStyle",lineSetting(i),"Color",colourSetting(i),"DisplayName",nameSetting(i))
    subplot(2,4,5)
    plot(t{i},prediction(2,:),"LineStyle",lineSetting(i),"Color",colourSetting(i),"DisplayName",nameSetting(i))
    subplot(2,4,3)
    plot(t{i},prediction(3,:),"LineStyle",lineSetting(i),"Color",colourSetting(i),"DisplayName",nameSetting(i))
    subplot(2,4,2)
    plot(t{i},prediction(4,:),"LineStyle",lineSetting(i),"Color",colourSetting(i),"DisplayName",nameSetting(i))
    subplot(2,4,6)
    plot(t{i},prediction(5,:),"LineStyle",lineSetting(i),"Color",colourSetting(i),"DisplayName",nameSetting(i))
    subplot(2,4,4)
    plot(t{i},prediction(6,:),"LineStyle",lineSetting(i),"Color",colourSetting(i),"DisplayName",nameSetting(i))
    subplot(2,4,7)
    plot(t{i},vmPrediction{i},"LineStyle",lineSetting(i),"Color",colourSetting(i),"DisplayName",nameSetting(i))
end
MeRE_MaRE = [vmMeRE_History, vmMaRE_History];
 excel = [MeRE_MaRE(1,:) MeRE_MaRE(2,:) MeRE_MaRE(3,:) ];
%% functions
function vmStress = calculateVonMisesStress(stress)
s11 = stress(1,:);
s22 = stress(2,:);
s33 = stress(3,:);
s12 = stress(4,:);
s23 = stress(5,:);
s13 = stress(6,:);
vmStress = sqrt(((s11-s22).^2 + (s22-s33).^2 + (s33-s11).^2 + 6.*(s12.^2+s23.^2+s13.^2))/2);
end

function [] = plottingInput(t, input, target, targetColour, targetLineType, targetMarker,fontSizeSetting,fontNameSetting)
figure
subplot(2,4,1)
% scatter(t(1:3:end),target(1,1:3:end),25,"Marker","x","MarkerEdgeColor","r","LineWidth",0.75,"DisplayName","Target")
plot(t(1:3:end),target(1,1:3:end),"Color",targetColour,"LineStyle",targetLineType,"Marker",targetMarker,"MarkerEdgeColor","auto","DisplayName","Target")
set(gca,'TickLabelInterpreter', 'latex','fontsize',fontSizeSetting,'FontName',fontNameSetting);
ylabel('$\sigma_{11}$ [MPa]','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
xlabel('t [-]','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
axis square
grid off
hold on

subplot(2,4,5)
% scatter(t(1:3:end),target(2,1:3:end),25,"Marker","x","MarkerEdgeColor","r","LineWidth",0.75,"DisplayName","Target")
plot(t(1:3:end),target(2,1:3:end),"Color",targetColour,"LineStyle",targetLineType,"Marker",targetMarker,"MarkerEdgeColor","auto","DisplayName","Target")
set(gca,'TickLabelInterpreter', 'latex','fontsize',fontSizeSetting,'FontName',fontNameSetting);
ylabel('$\sigma_{22}$ [MPa]','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
xlabel('t [-]','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
% title('$\sigma_{22}$ Response','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
axis square
grid off
hold on

subplot(2,4,3)
% scatter(t(1:3:end),target(3,1:3:end),25,"Marker","x","MarkerEdgeColor","r","LineWidth",0.75,"DisplayName","Target")
plot(t(1:3:end),target(3,1:3:end),"Color",targetColour,"LineStyle",targetLineType,"Marker",targetMarker,"MarkerEdgeColor","auto","DisplayName","Target")
set(gca,'TickLabelInterpreter', 'latex','fontsize',fontSizeSetting,'FontName',fontNameSetting);
ylabel('$\sigma_{33}$ [MPa]','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
xlabel('t [-]','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
% title('$\sigma_{33}$ Response','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
axis square
grid off
hold on

subplot(2,4,2)
% scatter(t(1:3:end),target(4,1:3:end),25,"Marker","x","MarkerEdgeColor","r","LineWidth",0.75,"DisplayName","Target")
plot(t(1:3:end),target(4,1:3:end),"Color",targetColour,"LineStyle",targetLineType,"Marker",targetMarker,"MarkerEdgeColor","auto","DisplayName","Target")
set(gca,'TickLabelInterpreter', 'latex','fontsize',fontSizeSetting,'FontName',fontNameSetting);
ylabel('$\sigma_{12}$ [MPa]','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
xlabel('t [-]','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
% title('$\sigma_{12}$ Response','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
axis square
grid off
hold on

subplot(2,4,6)
% scatter(t(1:3:end),target(5,1:3:end),25,"Marker","x","MarkerEdgeColor","r","LineWidth",0.75,"DisplayName","Target")
plot(t(1:3:end),target(5,1:3:end),"Color",targetColour,"LineStyle",targetLineType,"Marker",targetMarker,"MarkerEdgeColor","auto","DisplayName","Target")
set(gca,'TickLabelInterpreter', 'latex','fontsize',fontSizeSetting,'FontName',fontNameSetting);
ylabel('$\sigma_{23}$ [MPa]','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
xlabel('t [-]','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
% title('$\sigma_{23}$ Response','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
axis square
grid off
hold on

subplot(2,4,4)
% scatter(t(1:3:end),target(6,1:3:end),25,"Marker","x","MarkerEdgeColor","r","LineWidth",0.75,"DisplayName","Target")
plot(t(1:3:end),target(6,1:3:end),"Color",targetColour,"LineStyle",targetLineType,"Marker",targetMarker,"MarkerEdgeColor","auto","DisplayName","Target")
set(gca,'TickLabelInterpreter', 'latex','fontsize',fontSizeSetting,'FontName',fontNameSetting);
ylabel('$\sigma_{13}$ [MPa]','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
xlabel('t [-]','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting) 
% title('$\sigma_{13}$ Response','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
axis square
grid off
hold on

subplot(2,4,7)
vmTarget = calculateVonMisesStress(target);
% scatter(t(1:3:end),vmTarget(1,1:3:end),25,"Marker","x","MarkerEdgeColor","r","LineWidth",0.75,"DisplayName","Target")
plot(t(1:3:end),vmTarget(1,1:3:end),"Color",targetColour,"LineStyle",targetLineType,"Marker",targetMarker,"MarkerEdgeColor","auto","DisplayName","Target")
set(gca,'TickLabelInterpreter', 'latex','fontsize',fontSizeSetting,'FontName',fontNameSetting);
ylabel('$\sigma_{V}$ [MPa]','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
xlabel('t [-]','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
% title('Calculated Von Mises Stress','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
legend("Location","northeastoutside",'interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
axis square
grid off
hold on

subplot(2,4,8)
a11 = input(1,1); a22 = input(2,1); a33 = input(3,1);
a12 = input(4,1); a13 = input(5,1); a23 = input(6,1);
A = [a11, a12, a13;
    a12, a22, a23;
    a13, a23, a33];
disp(A)

[v, lam] = eig(A);

[X,Y,Z] = ellipsoid(0,0,0,lam(1,1),lam(2,2),lam(3,3),100);

ROT = v*[X(:) Y(:) Z(:)]';
X = reshape(ROT(1,:),size(X));
Y = reshape(ROT(2,:),size(Y));
Z = reshape(ROT(3,:),size(Z));

hold on

s1 = surf(X,Y,-ones(101),'FaceColor',[0 0 0]);
s2 = surf(Z+1,Y,zeros(101),'FaceColor',[0 0 0]);
s3 = surf(X,Z+1,zeros(101),'FaceColor',[0 0 0]);
rotate(s2,[0 -1 0], 90,[1,0,0])
rotate(s3,[1 0 0], 90,[0,1,0])
s4 = surf(X,Y,Z);

s1.EdgeColor = 'none';
s2.EdgeColor = 'none';
s3.EdgeColor = 'none';
s4.EdgeColor = 'none';

axis equal
grid off

xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
xticks([-1 -0.5 0 0.5 1])
yticks([-1 -0.5 0 0.5 1])
zticks([-1 -0.5 0 0.5 1])

vf = input(7,1);

title({append("Volume Fraction = ",num2str(vf)),"Orientation tensor"},'interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
xlabel('$e_{1}$','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
ylabel('$e_{2}$','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)
zlabel('$e_{3}$','interpreter','latex','fontsize',fontSizeSetting,'FontName',fontNameSetting)

view(-45,20)
camlight('headlight')
end

function [xRotated, yRotated, r] = preformPrincipalRotation(x,y,rotate2PrincipalAxis,nInputs,nDir)
if rotate2PrincipalAxis > 0
    timeSteps = size(x,2);
    vf = x(7,1);

    a11 = x(1,1); a22 = x(2,1);  a33 = x(3,1);
    a12 = x(4,1); a13 = x(5,1);  a23 = x(6,1);
    tmpA = [a11 a12 a13; a12 a22 a23; a13 a23 a33];
    [tmpR,~] = eig(tmpA);
    tmpR = tmpR';
    switch nDir
    case 1
        tmpR = tmpR([1 2 3],:);
    case 2
        tmpR = tmpR([1 3 2],:);
    case 3
        tmpR = tmpR([2 1 3],:);
    case 4
        tmpR = tmpR([2 3 1],:);
    case 5
        tmpR = tmpR([3 1 2],:);
    case 6
        tmpR = tmpR([3 2 1],:);
    end
    r = tmpR;

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

    if nInputs == 13
        xRotated = [repmat(prop, timeSteps, 1) strainStress(:,1:6)]';
    elseif nInputs == 9
        xRotated = [repmat(prop, timeSteps, 1) strainStress(:,1:6)]';
        xRotated = xRotated([2,3,7:13],:);
    end
    yRotated = strainStress(:,7:12)';
    r = r(:,:,1);
else
    xRotated = x;
    yRotated = y;
    r = [];
end
end

function yRotated = rotateBackToOrginalCoordinateFromPrincipal(y,r)
if isempty(r)
    yRotated = y;
else
    r = r';
    timeSteps = size(y,2);

    s11 = y(1,:); s22 = y(2,:); s33 = y(3,:);
    s12 = y(4,:); s23 = y(5,:); s13 = y(6,:);
    tmpStress = [s11; s12; s13; s12; s22; s23; s13; s23; s33];
    tmpStress = reshape(tmpStress,[3,3,timeSteps]);
    r = repmat(r,1,1,timeSteps);

    stressRotated = pagemtimes(pagemtimes(r,tmpStress),'none',r,'transpose');
    stressRotated = reshape(stressRotated,[9,size(stressRotated,3)])';

    yRotated = [stressRotated(:,1), stressRotated(:,5), stressRotated(:,9), stressRotated(:,2), stressRotated(:,6), stressRotated(:,3)]';   
end
end
