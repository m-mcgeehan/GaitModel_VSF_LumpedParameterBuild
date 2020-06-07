%Code to:
% 1) Read in static MoCap data from v3d pipeline
% 2) align v3d and simscape coordinate systems
% 3) Filter kinematics and calculate average static joint angles
% 4) Downsample and filter force data
% 5) Generate limb length estimates (subjx.seg) to scale model geometry in "ScaleBoneGeometry.m"
% 6) Create data structure with:
%     a) kinmeatics (subjx.km)
%     b) forces (subjx.forces)
%     c) framerate information (subjx.framerate)
%     d) segment length estimates (subjx.seg)
% 
% Data in: matlab variable generated from v3d pipeline
% Data out: subject data structure with the contents listed above
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
static.km.Rankle = RAnkleAngleFilt;
static.km.Lankle = LAnkleAngleFilt;
static.km.Rhip = RHipAngleFilt;
static.km.Lhip = LHipAngleFilt;
static.km.Rknee = RKneeAngleFilt;
static.km.Lknee = LKneeAngleFilt;
static.km.PelvisAngle = PelvisAngleFilt;
static.km.PelvisCOMPos = PelvisCOMPos;
static.forces.FP1 = FP1filt;
static.forces.FP2 = FP2filt;
static.framerate = FRAME_RATE;
static.trial = extractAfter(string(FILE_NAME{1}),"VSF"); 

% static.forces.FP1COP = downsample(FP1COP{1},5); %Need to zero COP outside of stance
% static.forces.FP2COP = downsample(FP2COP{1},5);

%% Generate the LCS offsets for each segment (align each segment CS with GCS)
static.GCSoffset.footL = mean(GCS_footL{1});
static.GCSoffset.footR = mean(GCS_footR{1});
static.GCSoffset.shankL = mean(GCS_shankL{1});
static.GCSoffset.shankR = mean(GCS_shankR{1});
static.GCSoffset.legL = mean(GCS_legL{1});
static.GCSoffset.legR = mean(GCS_legR{1});
static.GCSoffset.pelvis = mean(GCS_pelvis{1});

%% Generate limb length estimates
%Format marker variables
LASIavg = x*(mean(LASI{1}(:,1:3))');
RASIavg = x*(mean(RASI{1}(:,1:3))');
LGTRavg = x*(mean(LGTR{1}(:,1:3))');
RGTRavg = x*(mean(RGTR{1}(:,1:3))');
LHCEavg = x*(mean(LHCE{1}(:,1:3))');
RHCEavg = x*(mean(RHCE{1}(:,1:3))');
LHLAavg = x*(mean(LHLA{1}(:,1:3))');
RHLAavg = x*(mean(RHLA{1}(:,1:3))');
LHMEavg = x*(mean(LHME{1}(:,1:3))');
RHMEavg = x*(mean(RHME{1}(:,1:3))');
LLEPavg = x*(mean(LLEP{1}(:,1:3))');
RLEPavg = x*(mean(RLEP{1}(:,1:3))');
LMEPavg = x*(mean(LMEP{1}(:,1:3))');
RMEPavg = x*(mean(RMEP{1}(:,1:3))');
LPATavg = x*(mean(LPAT{1}(:,1:3))');
RPATavg = x*(mean(RPAT{1}(:,1:3))');
LPSIavg = x*(mean(LPSI{1}(:,1:3))');
RPSIavg = x*(mean(RPSI{1}(:,1:3))');
LTCEavg = x*(mean(LTCE{1}(:,1:3))');
RTCEavg = x*(mean(RTCE{1}(:,1:3))');
LTLAavg = x*(mean(LTLA{1}(:,1:3))');
RTLAavg = x*(mean(RTLA{1}(:,1:3))');
LTMEavg = x*(mean(LTME{1}(:,1:3))');
RTMEavg = x*(mean(RTME{1}(:,1:3))');

%Create virtual markers
AJCL = [.5*(LHLAavg(1)+LHMEavg(1)),.5*(LHLAavg(2)+LHMEavg(2)),.5*(LHLAavg(3)+LHMEavg(3))];
KJCL = [.5*(LLEPavg(1)+LMEPavg(1)),.5*(LLEPavg(2)+LMEPavg(2)),.5*(LLEPavg(3)+LMEPavg(3))];
KJCR = [.5*(RLEPavg(1)+RMEPavg(1)),.5*(RLEPavg(2)+RMEPavg(2)),.5*(RLEPavg(3)+RMEPavg(3))];
ASISdist = (abs(RASIavg(1) - LASIavg(1))) + (abs(RASIavg(2) - LASIavg(2))) + (abs(RASIavg(3) - LASIavg(3)));
PSISdist = (abs(RPSIavg(1) - LPSIavg(1))) + (abs(RPSIavg(2) - LPSIavg(2))) + (abs(RPSIavg(3) - LPSIavg(3)));
ASISmidpt = .5*(RASIavg + LASIavg);
PSISmidpt = .5*(RPSIavg + LPSIavg);
HJCL = [LASIavg - (.3*ASISdist), LASIavg + (.14*ASISdist), LASIavg - (.22*ASISdist)]; %Bell and Brand regression model
HJCR = [RASIavg - (.3*ASISdist), RASIavg - (.14*ASISdist), RASIavg - (.22*ASISdist)]; %Bell and Brand regression model
ground = [0,0,0];

%Estimate segment lengths
footXdist = abs(LTCEavg(1) - LHCEavg(1))*100; %Convert from m to cm (matches simscape variable input units)
footYdist = abs(LTLAavg(2) - LTMEavg(2))*100;
footZdist = abs(AJCL(3) - ground(3))*100;

legLengthL = abs(HJCL(3) - KJCL(3))*100;
legLengthR = abs(HJCR(3) - KJCR(3))*100;

shankLengthL = abs(KJCL(3) - AJCL(3))*100;

pelvisXdist = abs(ASISmidpt(1) - PSISmidpt(1))*100;
pelvisYdist = abs(LGTRavg(2)  - RGTRavg(2))*100;
pelvisZdist = (abs(PSISmidpt(3) + - HJCL(3)))*100; 
ASISYdist = (abs(RASIavg(2) - LASIavg(2)))*100;
PSISYdist = (abs(RPSIavg(2) - LPSIavg(2)))*100;

%Save semgent lengths to data structure
static.seg.footX = footXdist;
static.seg.footY = footYdist;
static.seg.footZ = footZdist;
static.seg.legLengthL = legLengthL;
static.seg.legLengthR = legLengthR;
static.seg.shankLengthL = shankLengthL;
static.seg.pelvisX = pelvisXdist;
static.seg.pelvisY = pelvisYdist;
static.seg.pelvisZ = pelvisZdist;
static.seg.ASISYdist = ASISYdist;
static.seg.PSISYdist = PSISYdist;

clearvars -except static

