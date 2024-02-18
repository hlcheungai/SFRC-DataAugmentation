clear; clc;

load("RVE_all_data.mat")
saveName = "RVE_all_data_RR_20_5";
nRandomRotation = 20;
skipSave = 0;

newX_train = cell(length(X_train)*nRandomRotation,1);
newY_train = cell(length(X_train)*nRandomRotation,1);
newX_valid = cell(length(X_valid)*nRandomRotation,1);
newY_valid = cell(length(X_valid)*nRandomRotation,1);
newX_test = cell(length(X_test)*nRandomRotation,1);
newY_test = cell(length(X_test)*nRandomRotation,1);


for j = 1:length(X_train)
    x = X_train{j};
    y = Y_train{j};
    for k = 1:nRandomRotation
        tmpR = createRandomRotationMatrix();
        [xRotated, yRotated] = preformRotation(x,y,tmpR);

        newX_train{length(X_train)*(k-1)+j} = xRotated;
        newY_train{length(X_train)*(k-1)+j} = yRotated;
    end
end

for j = 1:length(X_valid)
    x = X_valid{j};
    y = Y_valid{j};
    for k = 1:nRandomRotation
        tmpR = createRandomRotationMatrix();
        [xRotated, yRotated] = preformRotation(x,y,tmpR);

        newX_valid{length(X_valid)*(k-1)+j} = xRotated;
        newY_valid{length(X_valid)*(k-1)+j} = yRotated;
    end
end

for j = 1:length(X_test)
    x = X_test{j};
    y = Y_test{j};
    for k = 1:nRandomRotation
        tmpR = createRandomRotationMatrix();
        [xRotated, yRotated] = preformRotation(x,y,tmpR);

        newX_test{length(X_test)*(k-1)+j} = xRotated;
        newY_test{length(X_test)*(k-1)+j} = yRotated;
    end
end

X_train = newX_train;
Y_train = newY_train;
X_valid = newX_valid;
Y_valid = newY_valid;
X_test = newX_test;
Y_test = newY_test;

if skipSave == 0
    save(saveName,"X_train","Y_train","X_valid","Y_valid","X_test","Y_test")
end

function r = createRandomRotationMatrix()
% Random sampling of rotation parameters.
theta = 2*pi*rand();
phi = 2*pi*rand();
z = rand();

% Generate rotation matrix around z-axis.
R = [cos(theta) sin(theta) 0;
    -sin(theta) cos(theta) 0;
    0 0 1];

% Generate mirrored householder matrix.
v = [cos(phi)*sqrt(z); sin(phi)*sqrt(z); sqrt(1-z)];
P = 2*v*v' - eye(3);

% Total rotation matrix which is sampled uniformly from SO(3).
r = P*R;
end

function [xRotated, yRotated] = preformRotation(x,y,r)
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

strainRaotated = pagemtimes(pagemtimes(r,tmpStrain),'none',r,'transpose');
strainRaotated = reshape(strainRaotated,[9,size(strainRaotated,3)])';

stressRotated = pagemtimes(pagemtimes(r,tmpStress),'none',r,'transpose');
stressRotated = reshape(stressRotated,[9,size(stressRotated,3)])';

prop = [ARotated(1,1) ARotated(2,2) ARotated(3,3) ARotated(1,2) ARotated(1,3) ARotated(2,3) vf];

%strainStress = [strainRaotated(:,1), strainRaotated(:,5), strainRaotated(:,9), strainRaotated(:,2), strainRaotated(:,6), strainRaotated(:,3),...
%    stressRotated(:,1), stressRotated(:,5), stressRotated(:,9), stressRotated(:,2), stressRotated(:,6), stressRotated(:,3)];
strainStress = [strainRaotated(:,1), strainRaotated(:,5), strainRaotated(:,9), strainRaotated(:,2)*2, strainRaotated(:,6)*2, strainRaotated(:,3)*2,...
    stressRotated(:,1), stressRotated(:,5), stressRotated(:,9), stressRotated(:,2), stressRotated(:,6), stressRotated(:,3)];
xRotated = [repmat(prop, timeSteps, 1) strainStress(:,1:6)]';
yRotated = strainStress(:,7:12)';
vmOrginal =  calculateVonMisesStress(y);
vmRotated = calculateVonMisesStress(yRotated);
error = max(abs(vmOrginal-vmRotated))
if error > 10e-8
    disp("Check rotation matrix!")
    return
end
end

function vonMisesStress = calculateVonMisesStress(stress)
s11 = stress(1,:);
s22 = stress(2,:);
s33 = stress(3,:);
s12 = stress(4,:);
s23 = stress(5,:);
s13 = stress(6,:);
vonMisesStress = sqrt(((s11-s22).^2 + (s22-s33).^2 + (s33-s11).^2 + 6.*(s12.^2+s23.^2+s13.^2))/2);
end
