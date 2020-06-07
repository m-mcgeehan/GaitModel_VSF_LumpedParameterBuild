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

LAnkleAngleKinFoot = LAnkleAngleKinFoot{1}';
for i = 1:length(LAnkleAngleKinFoot(1,:))
    LAnkleAnglerot(i,:) = x*(LAnkleAngleKinFoot(:,i));
end
LAnkleAngleKinFoot = LAnkleAnglerot;
clear LAnkleAnglerot

LHipAngle = LHipAngle{1}';
for i = 1:length(LHipAngle(1,:))
    LHipAnglerot(i,:) = x*(LHipAngle(:,i));
end
LHipAngle = LHipAnglerot;
clear LHipAnglerot 

LKneeAngle = LKneeAngle{1}';
for i = 1:length(LKneeAngle(1,:))
    LKneeAnglerot(i,:) = x*(LKneeAngle(:,i));
end
LKneeAngle = LKneeAnglerot;
clear LKneeAnglerot

PelvisAngle = PelvisAngle{1}';
for i = 1:length(PelvisAngle(1,:))
    PelvisAnglerot(i,:) = x*(PelvisAngle(:,i));
end
PelvisAngle = PelvisAnglerot;
clear PelvisAnglerot

PelvisCOMPos = PelvisCOMPos{1}';
for i = 1:length(PelvisCOMPos(1,:))
    PelvisCOMPosrot(i,:) = x*(PelvisCOMPos(:,i));
end
PelvisCOMPos = PelvisCOMPosrot;
clear PelvisCOMPosrot
RAnkleAngle = RAnkleAngle{1}';
for i = 1:length(RAnkleAngle(1,:))
    RAnkleAnglerot(i,:) = x*(RAnkleAngle(:,i));
end
RAnkleAngle = RAnkleAnglerot;
clear RAnkleAnglerot

RAnkleAngleKinFoot = RAnkleAngleKinFoot{1}';
for i = 1:length(RAnkleAngleKinFoot(1,:))
    RAnkleAnglerot(i,:) = x*(RAnkleAngleKinFoot(:,i));
end
RAnkleAngleKinFoot = RAnkleAnglerot;
clear RAnkleAnglerot

RHipAngle = RHipAngle{1}';
for i = 1:length(RHipAngle(1,:))
    RHipAnglerot(i,:) = x*(RHipAngle(:,i));
end
RHipAngle = RHipAnglerot;
clear RHipAnglerot

RKneeAngle = RKneeAngle{1}';
for i = 1:length(RKneeAngle(1,:))
    RKneeAnglerot(i,:) = x*(RKneeAngle(:,i));
end
RKneeAngle = RKneeAnglerot;
clear RKneeAnglerot

%% Downsample forces and filter forces and COP
FP1ds = downsample(FP1{1},5);
FP2ds = downsample(FP2{1},5);
% FP1COPds = downsample(FP1COP{1},5);
% FP2COPds = downsample(FP2COP{1},5);

order = 2;
fs = 200;
fc = 50;%50 hz low pass cut off
Wn_fp = 2*fc/fs; 
[b,a] = butter(order,Wn_fp); 
FP1filt = filtfilt(b,a,FP1ds); %symmetric lowpass filter of normal force
FP2filt = filtfilt(b,a,FP2ds);
% FP1COP = filtfilt(b,a,FP1COPds);
% FP2COP = filtfilt(b,a,FP2COPds);
clear fc Wn_fp a b

%Zero COP outside of stance
%WRITE CODE FOR THIS

%% Filter kinematics
fc = 6;%6 hz low pass cut off
Wn_fp = 2*fc/fs; 
[b,a] = butter(order,Wn_fp);
RHipAngleFilt = filtfilt(b,a,RHipAngle);
LHipAngleFilt = filtfilt(b,a,LHipAngle);
RKneeAngleFilt = filtfilt(b,a,RKneeAngle);
LKneeAngleFilt = filtfilt(b,a,LKneeAngle);
RAnkleAngleFilt = filtfilt(b,a,RAnkleAngleKinFoot);
LAnkleAngleFilt = filtfilt(b,a,LAnkleAngleKinFoot);
PelvisAngleFilt = filtfilt(b,a,PelvisAngle);
PelvisCOMPosFilt = filtfilt(b,a,PelvisCOMPos);

%% set up structure
subj1.km.Rankle = RAnkleAngleFilt;
subj1.km.Lankle = LAnkleAngleFilt;
subj1.km.Rhip = RHipAngleFilt;
subj1.km.Lhip = LHipAngleFilt;
subj1.km.Rknee = RKneeAngleFilt;
subj1.km.Lknee = LKneeAngleFilt;
subj1.km.PelvisAngle = PelvisAngleFilt;
subj1.km.PelvisCOMPos = PelvisCOMPos;
subj1.forces.FP1 = FP1filt;
subj1.forces.FP2 = FP2filt;
subj1.framerate = FRAME_RATE;
subj1.seg = s1static.seg;
subj1.trial = extractAfter(string(FILE_NAME{1}),"VSF"); 
% subj1.forces.FP1COP = downsample(FP1COP{1},5); %Need to zero COP outside of stance
% subj1.forces.FP2COP = downsample(FP2COP{1},5);

