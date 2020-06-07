%Code to:
% 1) read in subject kinematic data from visual 3d pipeline
% 2) align visual 3d coordinate system with simscape
% 3) filter kinematic and force data and downsample forces
% 4) create data structure to drive model (usable in "startup_GaitModel_VSF.m")

%% Add path for dependencies
addpath(genpath('matlab scripts')); 

%% CONTROLS
staticNormalization_on = 0; %1 to normalize dynamic kinematics to static, 0 to bypass

%% Load mat file
[fName,fPath] = uigetfile('.mat');
load(fullfile(fPath,fName));
static = VSFabc_static;

%% DYNAMIC MOCAP PROCESSING
%Align coordinate systems (simscape wrt v3d: x = y, y = -x, z = z)
x = [0 -1 0;...
    1 0 0;...
    0 0 1];

%Rotate marker data and create data structure
% markerNames = {'RTME' 'RTLA' 'RTH4' 'RTH3' 'RTH2' 'RTH1' 'RTCE' 'RPSI' 'RMEP' 'RLEP' 'RHME' 'RHLA' 'RHCE' 'RGTR' 'RASI' 'LTME'...
%     'LTLA' 'LTH4' 'LTH3' 'LTCE' 'LSK4' 'LSK3' 'LSK2' 'LSK1' 'LPSI' 'LMEP' 'LLEP' 'LHME' 'LHLA' 'LHCE' 'LGTR' 'LASI'...
%     'RKJC' 'RHJC_reg' 'RHJC_fxnl' 'RAJC_fxnl'  'LHJC_reg' 'LHJC_fxnl' 'LAJC_fxnl' 'VSF_Rankle' 'CODA_origin'}; 


markerNames = {'RTME' 'RTLA' 'RTH4' 'RTH3' 'RTH2' 'RTH1' 'RTCE' 'RPSI' 'RMEP' 'RLEP' 'RHME' 'RHLA' 'RHCE' 'RGTR' 'RASI' 'LTME'...
    'LTLA' 'LTH4' 'LTH3' 'LTH2' 'LTH1' 'LTCE' 'LSK4' 'LSK3' 'LSK2' 'LSK1' 'LPSI' 'LMEP' 'LLEP' 'LHME' 'LHLA' 'LHCE' 'LGTR' 'LASI'...
    'RKJC' 'RHJC_reg' 'RHJC_fxnl' 'RAJC_fxnl' 'LKJC' 'LHJC_reg' 'LHJC_fxnl' 'LAJC_fxnl' 'VSF_Rankle' 'CODA_origin'};

order = 2;
fs = 200;
fc = 6;%6 hz low pass cut off
Wn_km = 2*fc/fs; 
[b,a] = butter(order,Wn_km); 
time = linspace(0,(length(LASI{1}(:,1))/fs),length(LASI{1}))';
for i = 1:length(markerNames)
    marker = markerNames{i};
    eval(['markerData =',marker,';']);
    markerDataIdx = markerData{1}(:,1:3)';
    markerDataFilt = filtfilt(b,a,markerDataIdx');
    markerDataIdx = markerDataFilt';
    for k = 1:length(markerDataIdx(1,:))
        markerDataRot(k,:) = (x*(markerDataIdx(:,k))*100);
    end
    markerSet.(marker) = markerDataRot;
    markerSetTS.(marker) = timeseries(markerDataRot,time);
end

clear 'RTME' 'RTLA' 'RTH4' 'RTH3' 'RTH2' 'RTH1' 'RTCE' 'RPSI' 'RMEP' 'RLEP' 'RHME' 'RHLA' 'RHCE' 'RGTR' 'RASI' 'LTME'...
    'LTLA' 'LTH4' 'LTH3' 'LTCE' 'LSK4' 'LSK3' 'LSK2' 'LSK1' 'LPSI' 'LMEP' 'LLEP' 'LHME' 'LHLA' 'LHCE' 'LGTR' 'LASI'...
    'RKJC' 'RHJC_reg' 'RHJC_fxnl' 'RAJC_fxnl' 'LKJC' 'LHJC_reg' 'LHJC_fxnl' 'LAJC_fxnl' 'VSF_Rankle' 'CODA_origin'
%% Create virtual markers and write them to the data structure
ASISmidpt = (markerSet.RASI+markerSet.LASI)/2;
PSISmidpt = (markerSet.RPSI+markerSet.LPSI)/2;
pelvisOrigin = ASISmidpt;
    ASISdist = norm(markerSet.RASI - markerSet.LASI);
HJCR = pelvisOrigin + (.1*[(11-(.063*static.legLength)) -1*(8+(.086*static.legLength)) (9-(.078*static.legLength))]);
HJCL = pelvisOrigin + (.1*[(11-(.063*static.legLength)) (8+(.086*static.legLength)) (9-(.078*static.legLength))]);
KJCR = (markerSet.RMEP + markerSet.RLEP)/2;
KJCL = (markerSet.LMEP + markerSet.LLEP)/2;
% FAJCL = markerSet.FAnkleL_landmark;
FAJCL = (markerSet.LHLA + markerSet.LHME)/2;

% HJCR = [markerSet.RASI(:,1) - (.36*ASISdist), markerSet.RASI(:,2) - (.14*ASISdist), markerSet.RASI(:,3) - (.22*ASISdist)];
% HJCL = [markerSet.LASI(:,1) - (.36*ASISdist), markerSet.LASI(:,2) + (.14*ASISdist), markerSet.LASI(:,3) - (.22*ASISdist)];
AJC_pros = [markerSet.RTCE(:,1) - markerSet.RHCE(:,1), markerSet.RHLA(:,2)+.4, markerSet.RHLA(:,3)+.62];
AJCL = (markerSet.LHLA + markerSet.LHME)/2;

markerSet.ASISmidpt = ASISmidpt;
markerSet.PSISmidpt = PSISmidpt;
markerSet.pelvisOrigin = pelvisOrigin;
% markerSet.HJCR = markerSet.FHipR_landmark;
% markerSet.HJCL = markerSet.FHipL_landmark;
markerSet.HJCR = HJCR;
markerSet.HJCL = HJCL;
% markerSet.KJCR = markerSet.FKneeR_landmark;
% markerSet.KJCL = markerSet.FKneeL_landmark;
markerSet.KJCR = KJCR;
markerSet.KJCL = KJCL;
% markerSet.AJCL = FAJCL;
markerSet.AJCL = AJCL;
markerSet.AJC_pros = AJC_pros;
markerSet.FAJC_pros = markerSet.LAJC_fxnl;

markerSetTS.ASISmidpt = timeseries(ASISmidpt,time);
markerSetTS.PSISmidpt = timeseries(PSISmidpt,time);
markerSetTS.pelvisOrigin = timeseries(pelvisOrigin,time);
markerSetTS.HJCR = timeseries(markerSet.HJCR,time);
markerSetTS.HJCL = timeseries(markerSet.HJCL,time);
markerSetTS.KJCR = timeseries(markerSet.KJCR,time);
markerSetTS.KJCL = timeseries(markerSet.KJCL,time);
markerSetTS.AJCL = timeseries(markerSet.AJCL,time);
markerSetTS.AJC_pros = timeseries(markerSet.AJC_pros,time);
markerSetTS.FACJ_pros = timeseries(markerSet.FAJC_pros,time);


clear i marker markerDataIdx markerData k markerDataRot 'FHipR_landmark' 'FHipL_landmark' 'FKneeL_landmark' 'FKneeR_landmark' 'FAnkleL_landmark' 'FAnkleR_landmark' 'LASI' 'RASI' 'LPSI'...
    'RPSI' 'LGTR' 'RGTR' 'LTH2' 'LTH3' 'LTH4' 'LPAT' 'RTH1' 'RTH2' 'RTH3' 'RTH4' 'RPAT' 'LMEP' 'LLEP' 'RMEP' 'RLEP' 'LSK1' 'LSK2' 'LSK3' 'LSK4'...
    'RSK1' 'RSK2' 'RSK3' 'RSK4' 'LHCE' 'RHCE' 'LHLA' 'RHLA' 'LHME' 'RHME' 'LTCE' 'RTCE' 'LTLA' 'RTLA' 'LTME' 'RTME'

%% v3d angles
pelvisAng_GCS_v3d = pelvisAng_GCS;
jointAngleNames = {'LAnkleAngle' 'LHipAngle' 'LKneeAngle' 'RAnkleAngle' 'RHipAngle' 'RKneeAngle' 'pelvisAng_GCS_v3d'};

order = 2;
fs = 200;
fc = 6;%6 hz low pass cut off
Wn_km = 2*fc/fs; 
[b,a] = butter(order,Wn_km); 
for i = 1:length(jointAngleNames)
    jointAngle = jointAngleNames{i};
    eval(['angle =',jointAngle,';']);
    angleDataIdx = angle{1}(:,1:3)';
    angleFilt = filtfilt(b,a,angleDataIdx');
    angleDataIdx = angleFilt';
    for k = 1:length(angleDataIdx(1,:))
        angleDataRot(k,:) = (x*(angleDataIdx(:,k)));
    end
    dynamic.km.v3d.(jointAngle) = angleDataRot;
    dynamic.km.v3d.TS.(jointAngle) = timeseries(angleDataRot,time);
end
% dynamic.forces.FP1COP = downsample(FP1COP{1},5); %Need to zero COP outside of stance
% dynamic.forces.FP2COP = downsample(FP2COP{1},5);

%% Estimate joint angle data and create data structure
numFrames = length(time);
%pelvis
for i = 1:numFrames
%     pelvis_j = (markerSet.LPSI(i,:) - markerSet.pelvisOrigin(i,:))/norm(markerSet.LPSI(i,:)-markerSet.pelvisOrigin(i,:));
%     pelvis_k = cross((markerSet.ASISmidpt(i,:) - markerSet.PSISmidpt(i,:))/(norm(markerSet.ASISmidpt(i,:)-markerSet.PSISmidpt(i,:))),pelvis_j);
    pelvis_j = (markerSet.LASI(i,:) - markerSet.pelvisOrigin(i,:))/norm((markerSet.LASI(i,:) - markerSet.pelvisOrigin(i,:)));
    pelvis_k = cross((markerSet.pelvisOrigin(i,:) - markerSet.PSISmidpt(i,:))/norm((markerSet.pelvisOrigin(i,:) - markerSet.PSISmidpt(i,:))),pelvis_j);
    pelvis_i = cross(pelvis_j,pelvis_k);
    pelvisR{i} = [pelvis_i', pelvis_j', pelvis_k'];
    if det(pelvisR{i}) < .95
    warning('determinant of pelvisR is below threshold. Determinant = ');
    display(det(pelvisR{i}));
    else
    end
end

%LegR
for i = 1:numFrames
    legR_k = (markerSet.HJCR(i,:) - markerSet.KJCR(i,:))/(norm(markerSet.HJCR(i,:) - markerSet.KJCR(i,:)));
    legR_i = cross(((markerSet.RMEP(i,:) - markerSet.RLEP(i,:))/norm(markerSet.RMEP(i,:) - markerSet.RLEP(i,:))),legR_k);
    legR_j = cross(legR_k,legR_i);
    legRR{i} = [legR_i', legR_j', legR_k'];
    if det(legRR{i}) < .95
    warning('determinant of legRR is below threshold. Determinant = ');
    display(det(legRR{i}));
    else
    end
end

%LegL
% for i = 1:numFrames
%     legL_k = (markerSet.HJCL(i,:) - markerSet.KJCL(i,:))/(norm(markerSet.HJCL(i,:) - markerSet.KJCL(i,:)));
%     legL_i = cross(((markerSet.LLEP(i,:) - markerSet.LMEP(i,:))/norm(markerSet.LLEP(i,:) - markerSet.LMEP(i,:))),legL_k);
%     legL_j = cross(legL_k,legL_i);
%     legLR{i} = [legL_i', legL_j', legL_k'];
%     if det(legLR{i}) < .95
%     warning('determinant of legLR is below threshold. Determinant = ');
%     display(det(legLR{i}));
%     else
%     end
% end
for i = 1:numFrames
    legL_k = (markerSet.LHJC_reg(i,:) - markerSet.LKJC(i,:))/(norm(markerSet.LHJC_reg(i,:) - markerSet.LKJC(i,:)));
    legL_i = cross(((markerSet.LLEP(i,:) - markerSet.LMEP(i,:))/norm(markerSet.LLEP(i,:) - markerSet.LMEP(i,:))),legL_k);
    legL_j = cross(legL_k,legL_i);
    legLR{i} = [legL_i', legL_j', legL_k'];
    if det(legLR{i}) < .95
    warning('determinant of legLR is below threshold. Determinant = ');
    display(det(legLR{i}));
    else
    end
end
%ShankR
for i = 1:numFrames
%     shankR_k = (markerSet.KJCR(i,:) - markerSet.FAJC_pros(i,:))/(norm(markerSet.KJCR(i,:) - markerSet.FAJC_pros(i,:)));
%     shankR_k = (markerSet.KJCR(i,:) - markerSet.AJC_pros(i,:))/(norm(markerSet.KJCR(i,:) - markerSet.AJC_pros(i,:)));
    shankR_k = (markerSet.KJCR(i,:) - markerSet.VSF_Rankle(i,:))/(norm(markerSet.KJCR(i,:) - markerSet.VSF_Rankle(i,:)));
    shankR_i = cross(((markerSet.RMEP(i,:) - markerSet.RLEP(i,:))/norm(markerSet.RMEP(i,:) - markerSet.RLEP(i,:))),shankR_k);
    shankR_j = cross(shankR_k,shankR_i);
    shankRR{i} = [shankR_i', shankR_j', shankR_k'];
    if det(shankRR{i}) < .95
    warning('determinant of shankRR is below threshold. Determinant = ');
    display(det(shankRR{i}));
    else
    end
end

%ShankL
for i = 1:numFrames
   shankL_k = (markerSet.LKJC(i,:) - markerSet.AJCL(i,:))/(norm(markerSet.LKJC(i,:) - markerSet.AJCL(i,:)));
%    shankL_k = (markerSet.LKJC(i,:) - markerSet.LAJC_fxnl(i,:))/(norm(markerSet.LKJC(i,:) - markerSet.LAJC_fxnl(i,:)));
   shankL_i = cross(((markerSet.LLEP(i,:) - markerSet.LMEP(i,:))/norm(markerSet.LLEP(i,:) - markerSet.LMEP(i,:))),shankL_k);
   shankL_j = cross(shankL_k,shankL_i);
   shankLR{i} = [shankL_i', shankL_j', shankL_k'];
   if det(shankLR{i}) < .95
    warning('determinant of shankLR is below threshold. Determinant = ');
    display(det(shankLR{i}));
   else
   end
end

%FootL 
for i = 1:numFrames
    footL_i = (markerSet.LTCE(i,:) - markerSet.LHCE(i,:))/norm(markerSet.LTCE(i,:) - markerSet.LHCE(i,:));
    footL_k = cross(footL_i,((markerSet.LTLA(i,:) - markerSet.LTME(i,:))/norm(markerSet.LTLA(i,:) - markerSet.LTME(i,:))));
    footL_j = cross(footL_k,footL_i);
    footLR{i} = [footL_i', footL_j', footL_k'];
    if det(footLR{i}) < .95
    warning('determinant of footRR is below threshold. Determinant = ');
    display(det(footLR{i}));
    else
    end
end

%Estimate segment angles wrt GCS (xyz euler sequence)
pelvisAng_GCS = zeros(length(time),3);
for i = 1:numFrames
    pelvisRi = pelvisR{i};
    pelvisAng_GCS(i,1:3) = wrapTo180(SpinCalc('DCMtoEA123',pelvisRi,.2,0));
end

legRAng_GCS = zeros(length(time),3);
for i = 1:numFrames
    legRRi = legRR{i};
    legRAng_GCS(i,1:3) = wrapTo180(SpinCalc('DCMtoEA123',legRRi,.2,0));
end

legLAng_GCS = zeros(length(time),3);
for i = 1:numFrames
    legLRi = legLR{i};
    legLAng_GCS(i,1:3) = wrapTo180(SpinCalc('DCMtoEA123',legLRi,.2,0));
end

shankRAng_GCS = zeros(length(time),3);
for i = 1:numFrames
    shankRRi = shankRR{i};
    shankRAng_GCS(i,1:3) = wrapTo180(SpinCalc('DCMtoEA123',shankRRi,.2,0));
end

shankLAng_GCS = zeros(length(time),3);
for i = 1:numFrames
    shankLRi = shankLR{i};
    shankLAng_GCS (i,1:3) = wrapTo180(SpinCalc('DCMtoEA123',shankLRi,.2,0));
end

footLAng_GCS = zeros(length(time),3);
for i = 1:numFrames
    footLRi = footLR{i};
    footLAng_GCS(i,1:3) = wrapTo180(SpinCalc('DCMtoEA123',footLRi,.2,0));
end

%estimate joint angles (distal wrt proximal)
pelvisRi = zeros(3,3);
legRRi = zeros(3,3);
RHipAngle = zeros (length(time),3);
for i = 1:numFrames
    pelvisRi = pelvisR{i};
    legRRi = legRR{i};
    RHipAngle(i,1:3) = wrapTo180(SpinCalc('DCMtoEA123',((inv(pelvisRi))*legRRi),.9,0));
end

legLRi = zeros(3,3);
LHipAngle = zeros (length(time),3);
for i = 1:numFrames
    pelvisRi = pelvisR{i};
    legLRi = legLR{i};
    LHipAngle(i,1:3) = wrapTo180(SpinCalc('DCMtoEA123',((inv(pelvisRi))*legLRi),.9,0));
end

shankRRi = zeros(3,3);
legRRi = zeros(3,3);
RKneeAngle = zeros (length(time),3);
for i = 1:numFrames
    legRRi = legRR{i};
    shankRRi = shankRR{i};
    RKneeAngle(i,1:3) = wrapTo180(SpinCalc('DCMtoEA123',((inv(legRRi))*shankRRi),.9,0));
end

shankLRi = zeros(3,3);
legLRi = zeros(3,3);
LKneeAngle = zeros (length(time),3);
for i = 1:numFrames
    legLRi = legLR{i};
    shankLRi = shankLR{i};
    LKneeAngle(i,1:3) = wrapTo180(SpinCalc('DCMtoEA123',((inv(legLRi))*shankLRi),.9,0));
end



% footRRi = zeros(3,3);
% shankRRi = zeros(3,3);
% RAnkleAngle = zeros (length(time),3);
% for i = 1:numFrames
%     shankRRi = shankRR{i};
%     footRRi = footRR{i};
%     RAnkleAngle(i,:) = wrapTo180(SpinCalc('DCMtoEA123',((inv(shankRRi))*footRRi),.9,0)); %precision tolerance violation 
% end

footLRi = zeros(3,3);
shankLRi = zeros(3,3);
LAnkleAngle = zeros(length(time),3);
for i = 1:numFrames
    shankLRi = shankLR{i};
    footLRi = footLR{i};
    LAnkleAngle(i,1:3) = wrapTo180(SpinCalc('DCMtoEA123',((inv(shankLRi))*footLRi),.9,0));
end

a = linspace(0,0,length(time));
RAnkleAngle = [a; a; a]';
%% Normalize joint angles to static motion capture (if user desires)
if staticNormalization_on == 1
    pelvisAng_GCS = [pelvisAng_GCS(:,1) (pelvisAng_GCS(:,2)-static.km.PelvisAngle(2)) pelvisAng_GCS(:,3)];
    RHipAngle = RHipAngle - static.km.Rhip;
    LHipAngle = LHipAngle - static.km.Lhip;
    RKneeAngle = RKneeAngle - static.km.Rknee;
    LKneeAngle = LKneeAngle - static.km.Lknee;
%     RKneeAngle = RKneeAngle - max(RKneeAngle);
%     LKneeAngle = LKneeAngle - max(LKneeAngle);
    
    LAnkleAngle = LAnkleAngle - static.km.Lankle;
end

dynamic.km.TS.pelvisAng_GCS = timeseries(pelvisAng_GCS,time);
dynamic.km.TS.RHipAngle = timeseries(RHipAngle,time);
dynamic.km.TS.LHipAngle = timeseries(LHipAngle,time);
dynamic.km.TS.RKneeAngle = timeseries(RKneeAngle,time);
dynamic.km.TS.LKneeAngle = timeseries(LKneeAngle,time);
dynamic.km.TS.LAnkleAngle = timeseries(LAnkleAngle,time);

%     
%% Scan through joint angle data and write to data structure
% jointAngleNames = {'LAnkleAngle' 'LHipAngle' 'LKneeAngle' 'pelvisAng_GCS' 'RHipAngle' 'RKneeAngle'};
% 
% for i = 1:length(jointAngleNames)
%     jointAngle = jointAngleNames{i};
%     eval(['JointAngleData =',jointAngle,';']);
%     jointAngles.(jointAngle) = JointAngleDataRot;
% end
% clear b a order fs fc Wn_km i jointAngle JointAngleData JointAngleDataIdx k JointAngleDataRot 'LAnkleAngle' 'LHipAngle' 'LKneeAngle' 'pelvisAng_GCS' 'RAnkleAngle' 'RHipAngle' 'RKneeAngle'

%% Downsample forces and filter forces and COP

order = 2;
fs = 1000;
fc = 40;%50 hz low pass cut off
Wn_fp = 2*fc/fs; 
[b,a] = butter(order,Wn_fp); 
FP1filt = filtfilt(b,a,FP1{1}); %symmetric lowpass filter of forces
FP2filt = filtfilt(b,a,FP2{1});

fc = 6;
Wn_fp = 2*fc/fs;
FP1COPfilt = filtfilt(b,a,FP1COP{1});
FP2COPfilt = filtfilt(b,a,FP2COP{1});

clear fc Wn_fp a b FP1 FP2 FP1COP FP2COP

FP1ds = downsample(FP1filt,10);
FP2ds = downsample(FP2filt,10);
FP1COPds = downsample(FP1COPfilt,10);
FP2COPds = downsample(FP2COPfilt,10);


%need to zero, filter, downsample, and rotate COP


for i =1:length(FP1ds(:,1))
    FP1(i,:) = x*(FP1ds(i,:)');
    FP2(i,:) = x*(FP2ds(i,:)');
    FP1COP(i,:) = x*(FP1COPds(i,:)');
    FP2COP(i,:) = x*(FP2COPds(i,:)'); 
end

clear FP1ds FP2ds FP1COPds FP2COPds


%% set up structure
dynamic.km.Lankle = LAnkleAngle;
dynamic.km.Rankle = RAnkleAngle;
dynamic.km.Rhip = RHipAngle;
dynamic.km.Lhip = LHipAngle;
dynamic.km.Rknee = RKneeAngle;
dynamic.km.Lknee = LKneeAngle;
dynamic.km.PelvisAngle = pelvisAng_GCS;
dynamic.km.PelvisCOMPos = markerSet.CODA_origin;
dynamic.markers = markerSet;
dynamic.markersTS = markerSetTS;
dynamic.forces.FP1 = FP1;
dynamic.forces.FP2 = FP2;
dynamic.forces.FP1COP = FP1COP;
dynamic.forces.FP2COP = FP2COP;
dynamic.framerate = FRAME_RATE;
dynamic.static = static;
dynamic.trial = fName; 

clearvars -except dynamic static fName

%% name and save data
trialInfo = replace(dynamic.trial,'-','_');
varName = convertStringsToChars(extractBefore(dynamic.trial,'.'));

eval([(varName) '=dynamic;']);

uisave((varName),(varName));
clear fName static dynamic varName trialInfo


