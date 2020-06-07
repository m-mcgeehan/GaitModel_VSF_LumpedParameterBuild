%Code to:
% 1) read in subject kinematic data from visual 3d pipeline
% 2) align visual 3d coordinate system with simscape
% 3) filter kinematic and force data and downsample forces
% 4) create data structure to drive model (usable in "startupLPMVSFWalker.m")

%% CONTROLS
staticNormalization_on = 1; %1 to normalize dynamic kinematics to static, 0 to bypass

%% DYNAMIC MOCAP PROCESSING
%Align coordinate systems (simscape wrt v3d: x = y, y = -x, z = z)
x = [0 1 0;...
    -1 0 0;...
    0 0 1];

%Rotate marker data and create data structure
markerNames = {'FHipR_landmark' 'FHipL_landmark' 'FKneeL_landmark' 'FKneeR_landmark' 'FAnkleL_landmark' 'FAnkleR_landmark' 'LASI' 'RASI' 'LPSI'...
    'RPSI' 'LGTR' 'RGTR' 'LTH2' 'LTH3' 'LTH4' 'LPAT' 'RTH1' 'RTH2' 'RTH3' 'RTH4' 'RPAT' 'LMEP' 'LLEP' 'RMEP' 'RLEP' 'LSK1' 'LSK2' 'LSK3' 'LSK4'...
    'RSK1' 'RSK2' 'RSK3' 'RSK4' 'LHCE' 'RHCE' 'LHLA' 'RHLA' 'LHME' 'RHME' 'LTCE' 'RTCE' 'LTLA' 'RTLA' 'LTME' 'RTME'};

order = 2;
fs = 200;
fc = 6;%50 hz low pass cut off
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
        markerDataRot(k,:) = x*(markerDataIdx(:,k));
    end
    markerSet.(marker) = markerDataRot;
    markerSetTS.(marker) = timeseries(markerDataRot,time);
end

clear i marker markerDataIdx markerData k markerDataRot 'FHipR_landmark' 'FHipL_landmark' 'FKneeL_landmark' 'FKneeR_landmark' 'FAnkleL_landmark' 'FAnkleR_landmark' 'LASI' 'RASI' 'LPSI'...
    'RPSI' 'LGTR' 'RGTR' 'LTH2' 'LTH3' 'LTH4' 'LPAT' 'RTH1' 'RTH2' 'RTH3' 'RTH4' 'RPAT' 'LMEP' 'LLEP' 'RMEP' 'RLEP' 'LSK1' 'LSK2' 'LSK3' 'LSK4'...
    'RSK1' 'RSK2' 'RSK3' 'RSK4' 'LHCE' 'RHCE' 'LHLA' 'RHLA' 'LHME' 'RHME' 'LTCE' 'RTCE' 'LTLA' 'RTLA' 'LTME' 'RTME'

%Rotate joint angle data and create data structure
jointAngleNames = {'LAnkleAngle' 'LHipAngle' 'LKneeAngle' 'pelvisAng_GCS' 'RAnkleAngle' 'RHipAngle' 'RKneeAngle'};

for i = 1:length(jointAngleNames)
    jointAngle = jointAngleNames{i};
    eval(['JointAngleData =',jointAngle,';']);
    JointAngleDataIdx = JointAngleData{1}(:,1:3)';
    JointAngleDataFilt = filtfilt(b,a,JointAngleDataIdx');
    JointAngleDataIdx = JointAngleDataFilt';
    for k = 1:length(JointAngleDataIdx(1,:))
        JointAngleDataRot(k,:) = x*(JointAngleDataIdx(:,k));
    end
    jointAngles.(jointAngle) = JointAngleDataRot;
end
clear b a order fs fc Wn_km i jointAngle JointAngleData JointAngleDataIdx k JointAngleDataRot 'LAnkleAngle' 'LHipAngle' 'LKneeAngle' 'pelvisAng_GCS' 'RAnkleAngle' 'RHipAngle' 'RKneeAngle'

%Pelvis COM position
PelvisCOMPos = PelvisCOMPos{1}';
for i = 1:length(PelvisCOMPos(1,:))
    PelvisCOMPosrot(i,:) = x*(PelvisCOMPos(:,i));
end
PelvisCOMPos = PelvisCOMPosrot;
clear PelvisCOMPosrot
%% Downsample forces and filter forces and COP

order = 2;
fs = 1000;
fc = 50;%50 hz low pass cut off
Wn_fp = 2*fc/fs; 
[b,a] = butter(order,Wn_fp); 
FP1filt = filtfilt(b,a,FP1{1}); %symmetric lowpass filter of forces
FP2filt = filtfilt(b,a,FP2{1});

clear fc Wn_fp a b

FP1ds = downsample(FP1filt,5);
FP2ds = downsample(FP2filt,5);

%% set up structure
subj1.km.Rankle = jointAngles.RAnkleAngle;
subj1.km.Lankle = jointAngles.LAnkleAngle;
subj1.km.Rhip = jointAngles.RHipAngle;
subj1.km.Lhip = jointAngles.LHipAngle;
subj1.km.Rknee = jointAngles.RKneeAngle;
subj1.km.Lknee = jointAngles.LKneeAngle;
subj1.km.PelvisAngle = jointAngles.pelvisAng_GCS;
subj1.km.PelvisCOMPos = PelvisCOMPos;
subj1.markers = markerSet;
subj1.markersTS = markerSetTS;
subj1.forces.FP1 = FP1ds;
subj1.forces.FP2 = FP2ds;
subj1.framerate = FRAME_RATE;
subj1.static = static;
subj1.trial = extractAfter(string(FILE_NAME{1}),"Visual3D\"); 



% subj1.forces.FP1COP = downsample(FP1COP{1},5); %Need to zero COP outside of stance
% subj1.forces.FP2COP = downsample(FP2COP{1},5);

clearvars -except subj1 static
