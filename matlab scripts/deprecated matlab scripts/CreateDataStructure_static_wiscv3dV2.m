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
%% Filter data
% Downsample forces and filter forces and COP
order = 2;
fs = 1000;
fc = 50;%50 hz low pass cut off
Wn_fp = 2*fc/fs; 
[b,a] = butter(order,Wn_fp); 
FP1filt = filtfilt(b,a,FP1{1}); %symmetric lowpass filter of normal force
FP2filt = filtfilt(b,a,FP2{1});
% FP1COP = filtfilt(b,a,FP1COPds);
% FP2COP = filtfilt(b,a,FP2COPds);
clear fc Wn_fp a b 

FP1ds = downsample(FP1filt,5);
FP2ds = downsample(FP2filt,5);
clear FP1filt FP2filt

%% Filter kinematics
fc = 6;%6 hz low pass cut off
Wn_km = 2*fc/fs; 
[b,a] = butter(order,Wn_km);
RHipAngleFilt = filtfilt(b,a,RHipAngle{1});
LHipAngleFilt = filtfilt(b,a,LHipAngle{1});
RKneeAngleFilt = filtfilt(b,a,RKneeAngle{1});
LKneeAngleFilt = filtfilt(b,a,LKneeAngle{1});
RAnkleAngleFilt = filtfilt(b,a,RAnkleAngle{1});
LAnkleAngleFilt = filtfilt(b,a,LAnkleAngle{1});
PelvisAngleFilt = filtfilt(b,a,pelvisAng_GCS{1});
PelvisCOMPosFilt = filtfilt(b,a,PelvisCOMPos{1});

%Segment orientations wrt GCS
LegLAng_GCS_filt = filtfilt(b,a,legLAng_GCS{1});
LegRAng_GCS_filt = filtfilt(b,a,legRAng_GCS{1});
ShankLAng_GCS_filt = filtfilt(b,a,shankLAng_GCS{1});
ShankRAng_GCS_filt = filtfilt(b,a,shankRAng_GCS{1});

clear fc fs Wn_km a b order RHipAngle LHipAngle RKneeAngle LKneeAngle RAnkleAngle LAnkleAngle pelvisAng_GCS PelvisCOMPos legLAng_GCS LegRAng_GCS legRAng_GCS shankLAng_GCS shankRAng_GCS

%% Align coordinate V3d and simscape coordinate systems
x = [0 1 0;...  %Align coordinate systems (simscape wrt v3d: x = y, y = -x, z = z)
    -1 0 0;...
    0 0 1];
%Joint angles   
LAnkleAngle = x*(mean(LAnkleAngleFilt)');
RAnkleAngle = x*(mean(RAnkleAngleFilt)');
LKneeAngle = x*(mean(LKneeAngleFilt)');
RKneeAngle = x*(mean(RKneeAngleFilt)');
LHipAngle = x*(mean(LHipAngleFilt)');
RHipAngle = x*(mean(RHipAngleFilt)');
PelvisAngle = x*(mean(PelvisAngleFilt)');
clear LAnkleAngleFilt RAnkleAngleFilt LKneeAngleFilt RKneeAngleFilt LHipAngleFilt RHipAngleFilt PelvisAngleFilt

%Pelvis COM pos
PelvisCOMPos = x*(mean(PelvisCOMPosFilt)'); clear PelvisCOMPosFilt

%Forces
FP1 = FP1{1}';
for i = 1:length(FP1(1,:))
    FP1rot(i,:) = x*(FP1(:,i));
end
FP1 = FP1rot;
clear FP1rot i

FP2 = FP2{1}';
for i = 1:length(FP2(1,:))
    FP2rot(i,:) = x*(FP2(:,i));
end
FP2 = FP2rot;
clear FP2rot i

FP1COP = FP1COP{1}';
for i = 1:length(FP1COP(1,:))
    FP1COProt(i,:) = x*(FP1COP(:,i));
end
FP1COP = FP1COProt*100;
clear FP1COProt i

FP2COP = FP2COP{1}';
for i = 1:length(FP2COP(1,:))
    FP2COProt(i,:) = x*(FP2COP(:,i));
end
FP2COP = FP2COProt*100;
clear FP2COProt i

%Markers
LASIavg = (x*(mean(LASI{1}(:,1:3))')*100);
RASIavg = (x*(mean(RASI{1}(:,1:3))')*100);
LPSIavg = (x*(mean(LPSI{1}(:,1:3))')*100);
RPSIavg = (x*(mean(RPSI{1}(:,1:3))')*100);
LGTRavg = (x*(mean(LGTR{1}(:,1:3))')*100);
RGTRavg = (x*(mean(RGTR{1}(:,1:3))')*100);
LTH1avg = (x*(mean(LTH1{1}(:,1:3))')*100);
LTH2avg = (x*(mean(LTH2{1}(:,1:3))')*100);
LTH3avg = (x*(mean(LTH3{1}(:,1:3))')*100);
LTH4avg = (x*(mean(LTH4{1}(:,1:3))')*100);
RTH1avg = (x*(mean(RTH1{1}(:,1:3))')*100);
RTH2avg = (x*(mean(RTH2{1}(:,1:3))')*100);
RTH3avg = (x*(mean(RTH3{1}(:,1:3))')*100);
RTH4avg = (x*(mean(RTH4{1}(:,1:3))')*100);
LLEPavg = (x*(mean(LLEP{1}(:,1:3))')*100);
RLEPavg = (x*(mean(RLEP{1}(:,1:3))')*100);
LMEPavg = (x*(mean(LMEP{1}(:,1:3))')*100);
RMEPavg = (x*(mean(RMEP{1}(:,1:3))')*100);
LPATavg = (x*(mean(LPAT{1}(:,1:3))')*100);
RPATavg = (x*(mean(RPAT{1}(:,1:3))')*100);
LSK1avg = (x*(mean(LSK1{1}(:,1:3))')*100);
LSK2avg = (x*(mean(LSK2{1}(:,1:3))')*100);
LSK3avg = (x*(mean(LSK3{1}(:,1:3))')*100);
LSK4avg = (x*(mean(LSK4{1}(:,1:3))')*100);
RSK1avg = (x*(mean(RSK1{1}(:,1:3))')*100);
RSK2avg = (x*(mean(RSK2{1}(:,1:3))')*100);
RSK3avg = (x*(mean(RSK3{1}(:,1:3))')*100);
RSK4avg = (x*(mean(RSK4{1}(:,1:3))')*100);
LHCEavg = (x*(mean(LHCE{1}(:,1:3))')*100);
RHCEavg = (x*(mean(RHCE{1}(:,1:3))')*100);
LHLAavg = (x*(mean(LHLA{1}(:,1:3))')*100);
RHLAavg = (x*(mean(RHLA{1}(:,1:3))')*100);
LHMEavg = (x*(mean(LHME{1}(:,1:3))')*100);
RHMEavg = (x*(mean(RHME{1}(:,1:3))')*100);
LTCEavg = (x*(mean(LTCE{1}(:,1:3))')*100);
RTCEavg = (x*(mean(RTCE{1}(:,1:3))')*100);
LTLAavg = (x*(mean(LTLA{1}(:,1:3))')*100);
RTLAavg = (x*(mean(RTLA{1}(:,1:3))')*100);
LTMEavg = (x*(mean(LTME{1}(:,1:3))')*100);
RTMEavg = (x*(mean(RTME{1}(:,1:3))')*100);

static.markers.LASI = LASIavg;
static.markers.RASI = RASIavg;
static.markers.LPSI = LPSIavg;
static.markers.RPSI = RPSIavg;
static.markers.LGTR = LGTRavg;
static.markers.RGTR = RGTRavg;
static.markers.LTH1 = LTH1avg;
static.markers.LTH2 = LTH2avg;
static.markers.LTH3 = LTH3avg;
static.markers.LTH4 = LTH4avg;
static.markers.RTH1 = RTH1avg;
static.markers.RTH2 = RTH2avg;
static.markers.RTH3 = RTH3avg;
static.markers.RTH4 = RTH4avg;
static.markers.LLEP = LLEPavg;
static.markers.RLEP = RLEPavg;
static.markers.LMEP = LMEPavg;
static.markers.RMEP = RMEPavg;
static.markers.LPAT = LPATavg;
static.markers.RPAT = RPATavg;
static.markers.LSK1 = LSK1avg;
static.markers.LSK2 = LSK2avg;
static.markers.LSK3 = LSK3avg;
static.markers.LSK4 = LSK4avg;
static.markers.RSK1 = RSK1avg;
static.markers.RSK2 = RSK2avg;
static.markers.RSK3 = RSK3avg;
static.markers.RSK4 = RSK4avg;
static.markers.LHCE = LHCEavg;
static.markers.RHCE = RHCEavg;
static.markers.LHLA = LHLAavg;
static.markers.RHLA = RHLAavg;
static.markers.LHME = LHMEavg;
static.markers.RHME = RHMEavg;
static.markers.LTCE = LTCEavg;
static.markers.RTCE = RTCEavg;
static.markers.LTLA = LTLAavg;
static.markers.RTLA = RTLAavg;
static.markers.LTME = LTMEavg;
static.markers.RTME = RTMEavg;

%Should filter marker coordinates at some point
clear LASI RASI LHCE RHCE LHLA RHLA LHME RHME LLEP RLEP LMEP RMEP LPSI RPSI LTCE RTCE LTLA RTLA LTME RTME LSK1 LSK2 LSK3 LSK4 LASI_processed

%Format virtual marks for joint centers
FAJCLPos_ShankL = (x*(mean(FAJCLPos_ShankL{1}(:,1:3))')*100);
FAJCRPos_ShankR = (x*(mean(FAJCRPos_ShankR{1}(:,1:3))')*100);
FHipLPos_LegL = (x*(mean(FHipLPos_LegL{1}(:,1:3))')*100);
FHipLPos_pelvis = (x*(mean(FHipLPos_pelvis{1}(:,1:3))')*100);
FHipRPos_LegR = (x*(mean(FHipRPos_LegR{1}(:,1:3))')*100);
FHipRPos_pelvis = (x*(mean(FHipRPos_pelvis{1}(:,1:3))')*100);
FKneeLPos_LegL = (x*(mean(FKneeLPos_LegL{1}(:,1:3))')*100);
FKneeLPos_ShankL = (x*(mean(FKneeLPos_ShankL{1}(:,1:3))')*100);
FKneeRPos_LegR = (x*(mean(FKneeRPos_LegR{1}(:,1:3))')*100);
FKneeRPos_ShankR = (x*(mean(FKneeRPos_ShankR{1}(:,1:3))')*100);
FootLPos_ShankL = (x*(mean(FootLPos_ShankL{1}(:,1:3))')*100);
FootRPos_ShankR = (x*(mean(FootRPos_ShankR{1}(:,1:3))')*100);
AJC_pros = [(RTCEavg(1)-RHCEavg(1)), RHLAavg(2)+.04, RHLAavg(3)+.0062]; 

%Create virtual markers
ASISdist = (abs(RASIavg(1) - LASIavg(1))) + (abs(RASIavg(2) - LASIavg(2))) + (abs(RASIavg(3) - LASIavg(3)));
PSISdist = (abs(RPSIavg(1) - LPSIavg(1))) + (abs(RPSIavg(2) - LPSIavg(2))) + (abs(RPSIavg(3) - LPSIavg(3)));
ASISmidpt = .5*(RASIavg + LASIavg);
PSISmidpt = .5*(RPSIavg + LPSIavg);
pelvisOrigin = (ASISmidpt+PSISmidpt)/2;
LSKmid = ((((LSK2avg+LSK3avg)/2) + (LSK4avg))/2); 
ground = [0,0,0];

%% Generate the LCS translational offsets for each segment and save to data structure
static.offsets.FAJCLPos_ShankL = FAJCLPos_ShankL;
static.offsets.FAJCRPos_ShankR = FAJCRPos_ShankR;
static.offsets.FAJCLPos_GCS = x*(mean(FAnkleL_landmark{1})'*100); 
static.offsets.FAJCRPos_GCS = x*(mean(FAnkleR_landmark{1})'*100); clear FAJCLPos_ShankL FAJCRPos_ShankR FAnkleL_landmark FAnkleR_landmark

static.offsets.FHipLPos_LegL = FHipLPos_LegL;
static.offsets.FHipRPos_LegR = FHipRPos_LegR;
static.offsets.FHipLPos_GCS = x*(mean(FHipL_landmark{1})'*100); 
static.offsets.FHipRPos_GCS = x*(mean(FHipR_landmark{1})'*100); clear FHipLPos_LegL FHipLPos_LegR FHipL_landmark FHipR_landmark

static.offsets.FKneeLPos_LegL = FKneeLPos_LegL;
static.offsets.FKneeRPos_LegR = FKneeRPos_LegR;
static.offsets.FKneeLPos_GCS = x*(mean(FKneeL_landmark{1})'*100); 
static.offsets.FKneeRPos_GCS = x*(mean(FKneeR_landmark{1})'*100); clear FKneeLPos_LegL FKneeRPos_LegR F FKneeL_landmark FKneeR_landmark

static.offsets.LegLDistalPos_GCS = x*(mean(LegLDistalPos{1})'*100);
static.offsets.LegLProxPos_GCS = x*(mean(LegLProxPos{1})'*100);
static.offsets.LegRDistalPos_GCS = x*(mean(LegRDistalPos{1})'*100);
static.offsets.LegRProxPos_GCS = x*(mean(LegRProxPos{1})'*100); 
static.offsets.LegLProxPos_pelvis = x*(mean(LegLPos_pelvis{1})'*100);
static.offsets.LegRProxPos_pelvis = x*(mean(LegRPos_pelvis{1})'*100); clear LegLDistalPos LegLProxPos LegRDistalPos  LegRProxPos LegLPos_pelvis LegRPos_pelvis

static.offsets.ShankLDistalPos_GCS = x*(mean(ShankLDistalPos{1})'*100);
static.offsets.ShankLProxPos_GCS = x*(mean(ShankLProxPos{1})'*100);
static.offsets.ShankRDistalPos_GCS = x*(mean(ShankRDistalPos{1})'*100);
static.offsets.ShankRProxPos_GCS = x*(mean(ShankRProxPos{1})'*100); 
static.offsets.ShankLProxPos_LegL = x*(mean(ShankLPos_LegL{1})'*100);
static.offsets.ShankRProxPos_LegR = x*(mean(ShankRPos_LegR{1})'*100); clear ShankLDistalPos ShankLProxPos ShankRDistalPos  ShankRProxPos ShankLPos_LegL ShankRPos_LegR

static.offsets.PelvisDistPos_GCS = x*(mean(PelivsRDistalPos{1})'*100);
static.offsets.PelvisProxPos_GCS = x*(mean(PelivsRProxPos{1})'*100); clear PelivsRDistalPos PelivsRProxPos

static.markers.FAJCL = static.offsets.FAJCLPos_GCS;
static.markers.FAJCR = static.offsets.FAJCLPos_GCS;
static.markers.FKJCL = static.offsets.FKneeLPos_GCS;
static.markers.FKJCR = static.offsets.FKneeRPos_GCS;
static.markers.FHJCL = static.offsets.FHipLPos_GCS;
static.markers.FHJCR = static.offsets.FHipRPos_GCS;

%% Generate the LCS rotational offsets for each segment and save to data structure
static.GCSOrientations.LegLAng_GCS = x*(mean(LegLAng_GCS_filt)');
static.GCSOrientations.LegRAng_GCS = x*(mean(LegRAng_GCS_filt)');
static.GCSOrientations.ShankLAng_GCS = x*(mean(ShankLAng_GCS_filt)');
static.GCSOrientations.ShankRAng_GCS = x*(mean(ShankRAng_GCS_filt)');
static.GCSOrientations.Pelvis_GCS = PelvisAngle; clear LegLAng_GCS_filt LegRAng_GCS_filt ShankLAng_GCS_filt ShankRAng_GCS_filt

%% Generate limb length estimates
%Estimate segment lengths
legLengthL = norm(static.offsets.FHipLPos_GCS - static.offsets.FKneeLPos_GCS); 
legLengthR = norm(static.offsets.FHipRPos_GCS - static.offsets.FKneeRPos_GCS); 

shankLengthL = norm(static.offsets.ShankLProxPos_GCS - static.offsets.FAJCLPos_GCS); %YOU ARE HERE
shankLengthR = norm(static.offsets.ShankRProxPos_GCS(3) - (LSKmid(3)));

pelvisXdist = norm(ASISmidpt - PSISmidpt);
pelvisYdist = norm(RASIavg - LASIavg);
pelvisZdist = (abs(PSISmidpt(3) - static.offsets.FHipLPos_GCS(3))); 
ASISYdist = (abs(RASIavg(2) - LASIavg(2)));
PSISYdist = (abs(RPSIavg(2) - LPSIavg(2)));

footLdistX = abs(LTCEavg(1) - LHCEavg(1)); %Convert from m to cm (matches simscape variable input units)
footLdistY = abs(LTLAavg(2) - LTMEavg(2));
footLdistZ = static.offsets.FAJCLPos_GCS(3) - ground(3);

%% set up data structure
static.km.Rankle = RAnkleAngle;
static.km.Lankle = LAnkleAngle;
static.km.Rhip = RHipAngle;
static.km.Lhip = LHipAngle;
static.km.Rknee = RKneeAngle;
static.km.Lknee = LKneeAngle;
static.km.PelvisAngle = PelvisAngle;
static.km.PelvisCOMPos = PelvisCOMPos;
static.forces.FP1 = FP1ds;
static.forces.FP2 = FP2ds;
static.framerate = FRAME_RATE;
static.trial = extractAfter(string(FILE_NAME{1}),"VSF\"); 
static.mass = mass{1};
static.height = height{1};
%% Save segment length estimates to data structure 
static.segLength.footLX = footLdistX;
static.segLength.footLY = footLdistY;
static.segLength.footLZ = footLdistZ;
static.segLength.legLengthL = legLengthL;
static.segLength.legLengthR = legLengthR;
static.segLength.shankLengthL = shankLengthL;
static.segLength.shankLengthR = shankLengthR;
static.segLength.pelvisX = pelvisXdist;
static.segLength.pelvisY = pelvisYdist;
static.segLength.pelvisZ = pelvisZdist;
static.segLength.ASISYdist = ASISYdist;
static.segLength.PSISYdist = PSISYdist;

clearvars -except static
