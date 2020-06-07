%Code to:
% 1) Read in static MoCap data from v3d pipeline
% 2) align v3d and simscape coordinate systems
% 3) Filter kinematics and calculate average static joint angles
% 4) Downsample and filter force data
% 5) Generate limb length estimates scale model geometry in "ScaleBoneGeometry.m"
% 6) Create data structure with:
%     a) Marker coordinates (static.markers)
    % b) segment angles wrt the GCS (static.GCSOrientations)
    % c) joint angles distal wrt proximal) (static.km)
    % d) joint center locations and limb translational offsets (static.offsets)
    % e) Generate limb length estimates (static.segLength) to scale model geometry in "ScaleBoneGeometry.m"
    % f) forces (static.forces)
    % g) frame rate information (static.framerate)
    % h) trial information (static.trial)
    % i) subject mass (static.mass)
    % j) subject height (static.height)
% 
% Data in: matlab variable generated from v3d pipeline
% Data out: data structure ('static') with the contents listed above
%% Add file paths for dependencies
addpath(genpath('matlab scripts')); 

%% Rotate markers and setup marker structure
x = [0 1 0;...  %rotation matrix to align coordinate systems (simscape wrt v3d: x = y, y = -x, z = z)
    -1 0 0;...
    0 0 1];

%Markers
LASIavg = (x*(mean(LASI{1}(:,1:3))')*100);
RASIavg = (x*(mean(RASI{1}(:,1:3))')*100);
LPSIavg = (x*(mean(LPSI{1}(:,1:3))')*100);
RPSIavg = (x*(mean(RPSI{1}(:,1:3))')*100);
LGTRavg = (x*(mean(LGTR{1}(:,1:3))')*100);
RGTRavg = (x*(mean(RGTR{1}(:,1:3))')*100);
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
% LPATavg = (x*(mean(LPAT{1}(:,1:3))')*100);
% RPATavg = (x*(mean(RPAT{1}(:,1:3))')*100);
LSK1avg = (x*(mean(LSK1{1}(:,1:3))')*100);
LSK2avg = (x*(mean(LSK2{1}(:,1:3))')*100);
LSK3avg = (x*(mean(LSK3{1}(:,1:3))')*100);
LSK4avg = (x*(mean(LSK4{1}(:,1:3))')*100);
% RSK1avg = (x*(mean(RSK1{1}(:,1:3))')*100);
% RSK2avg = (x*(mean(RSK2{1}(:,1:3))')*100);
% RSK3avg = (x*(mean(RSK3{1}(:,1:3))')*100);
% RSK4avg = (x*(mean(RSK4{1}(:,1:3))')*100);
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
LSKmid = ((((LSK2avg+LSK3avg)/2) + (LSK4avg))/2); 

static.markers.LASI = LASIavg;
static.markers.RASI = RASIavg;
static.markers.LPSI = LPSIavg;
static.markers.RPSI = RPSIavg;
static.markers.LGTR = LGTRavg;
static.markers.RGTR = RGTRavg;
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
% static.markers.LPAT = LPATavg;
% static.markers.RPAT = RPATavg;
static.markers.LSK1 = LSK1avg;
static.markers.LSK2 = LSK2avg;
static.markers.LSK3 = LSK3avg;
static.markers.LSK4 = LSK4avg;
static.markers.LSKmid = LSKmid;
% static.markers.RSK1 = RSK1avg;
% static.markers.RSK2 = RSK2avg;
% static.markers.RSK3 = RSK3avg;
% static.markers.RSK4 = RSK4avg;
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
static.markers.ground = [0,0,0];

clear LASI RASI LPSI RPSI LGTR RGTR LTH1 LTH2 LTH3 LTH4 RTH1 RTH2 RTH3 RTH4 LHCE RHCE LHLA RHLA LHME RHME LLEP RLEP...
    LMEP RMEP LPSI RPSI LTCE RTCE LTLA RTLA LTME RTME LSK1 LSK2 LSK3 LSK4 RSK1 RSK2 RSK3 RSK4 LASI_processed
%%
%format virtual markers
ASISmidpt = (static.markers.RASI+static.markers.LASI)/2;
PSISmidpt = (static.markers.RPSI+static.markers.LPSI)/2;
% pelvisOrigin = PSISmidpt;
pelvisOrigin = ASISmidpt;
    ASISdist = norm(static.markers.RASI - static.markers.LASI);
% HJCR = (x*(mean(FHipR_landmark{1}(:,1:3))')*100);
% HJCL = (x*(mean(FHipL_landmark{1}(:,1:3))')*100);
legLength = 10*(norm(ASISmidpt - static.markers.LMEP) + norm(static.markers.LMEP - static.markers.LHME));
HJCR = pelvisOrigin + (.1*[(11-(.063*legLength)); -1*(8+(.086*legLength)); (9-(.078*legLength))]);
HJCL = pelvisOrigin + (.1*[(11-(.063*legLength)); (8+(.086*legLength)); (9-(.078*legLength))]);
KJCR = (static.markers.RMEP + static.markers.RLEP)/2;
KJCL = (static.markers.LMEP + static.markers.LLEP)/2;
% KJCR = (x*(mean(FKneeR_landmark{1}(:,1:3))')*100);
% KJCL = (x*(mean(FKneeL_landmark{1}(:,1:3))')*100);
AJC_pros = [static.markers.RTCE(1) - static.markers.RHCE(1), static.markers.RHLA(2)+.4, static.markers.RHLA(3)+.62]';
AJCL = (static.markers.LHLA + static.markers.LHME)/2;

static.markers.ASISmidpt = ASISmidpt;
static.markers.PSISmidpt = PSISmidpt;
static.markers.pelvisOrigin = pelvisOrigin;
static.markers.HJCR = HJCR;
static.markers.HJCL = HJCL;
static.markers.KJCR = KJCR;
static.markers.KJCL = KJCL;
static.markers.AJC_pros = AJC_pros;
static.markers.AJCR = (x*(mean(FAnkleR_landmark{1}(:,1:3))')*100);
% static.markers.AJCL = (x*(mean(FAnkleL_landmark{1}(:,1:3))')*100);
static.markers.AJCL = AJCL;
clear ASIS midpt PSISmidpt pelvisOrigin HJCR HJCL KJCR KJCL AJC_pros AJCR AJCL 
%% Estimate kinematics
%pelvis
% pelvis_j = (static.markers.LPSI-static.markers.pelvisOrigin)/norm(static.markers.LPSI - static.markers.pelvisOrigin);
% pelvis_k = cross((static.markers.ASISmidpt-static.markers.PSISmidpt)/(norm(static.markers.ASISmidpt-static.markers.PSISmidpt)),pelvis_j);
pelvis_j = (static.markers.LASI - static.markers.pelvisOrigin)/norm(static.markers.LASI - static.markers.pelvisOrigin);
pelvis_k = cross((static.markers.pelvisOrigin-static.markers.PSISmidpt)/(norm(static.markers.pelvisOrigin-static.markers.PSISmidpt)),pelvis_j);
pelvis_i = cross(pelvis_j,pelvis_k);
pelvisR = [pelvis_i, pelvis_j, pelvis_k];

if det(pelvisR) < .95
warning('determinant of pelvisR is below threshold. Determinant = ');
display(det(pelvisR));
else
end

%Leg
legR_k = (static.markers.HJCR - static.markers.KJCR)/(norm(static.markers.HJCR - static.markers.KJCR));
legR_i = cross(((static.markers.RMEP - static.markers.RLEP)/norm(static.markers.RMEP - static.markers.RLEP)),legR_k);
legR_j = cross(legR_k,legR_i);
legRR = [legR_i, legR_j, legR_k]';

if det(legRR) < .95
warning('determinant of legRR is below threshold. Determinant = ');
display(det(legRR));
else
end

legL_k = (static.markers.HJCL - static.markers.KJCL)/(norm(static.markers.HJCL - static.markers.KJCL));
legL_i = cross(((static.markers.LLEP - static.markers.LMEP)/norm(static.markers.LLEP - static.markers.LMEP)),legL_k); 
legL_j = cross(legL_k,legL_i);
legLR = [legL_i, legL_j, legL_k]';
legLRdet = det(legLR);

if det(legLR) < .95
warning('determinant of legLR is below threshold. Determinant = ');
display(det(legLR));
else
end

%Shank
shankR_k = (static.markers.KJCR - static.markers.AJCR)/(norm(static.markers.KJCR - static.markers.AJCR));
shankR_i = cross(((static.markers.RMEP - static.markers.RLEP)/norm(static.markers.RMEP - static.markers.RLEP)),shankR_k);
shankR_j = cross(shankR_k,shankR_i);
shankRR = [shankR_i, shankR_j, shankR_k]';

if det(shankRR) < .95
warning('determinant of shankRR is below threshold. Determinant = ');
display(det(shankRR));
else
end

shankL_k = (static.markers.KJCL - static.markers.AJCL)/(norm(static.markers.KJCL - static.markers.AJCL));
shankL_i = cross(((static.markers.LLEP - static.markers.LMEP)/norm(static.markers.LLEP - static.markers.LMEP)),shankL_k);
shankL_j = cross(shankL_k,shankL_i);
shankLR = [shankL_i, shankL_j, shankL_k];

if det(shankLR) < .95
warning('determinant of shankLR is below threshold. Determinant = ');
display(det(shankLR));
else
end

%Foot
footR_i = (static.markers.RTCE - static.markers.RHCE)/norm(static.markers.RTCE - static.markers.RHCE);
footR_k = cross(footR_i,((static.markers.RTME - static.markers.RTLA)/norm(static.markers.RTME - static.markers.RTLA)));
footR_j = cross(footR_k,footR_i);
footRR = [footR_i, footR_j, footR_k]';

if det(footRR) < .95
warning('determinant of footRR is below threshold. Determinant = ');
display(det(footRR));
else
end

footL_i = (static.markers.LTCE - static.markers.LHCE)/norm(static.markers.LTCE - static.markers.LHCE);
footL_k = cross(footL_i,((static.markers.LTLA - static.markers.LTME)/norm(static.markers.LTLA - static.markers.LTME)));
footL_j = cross(footL_k,footL_i);
footLR = [footL_i, footL_j, footL_k]';

if det(footLR) < .95
warning('determinant of footLR is below threshold...Determinant =  ');
display(det(footLR));
else
end

% Calculate segment angles wrt GCS
pelvisAng_GCS = wrapTo180(SpinCalc('DCMtoEA123',pelvisR,.2,1));
legRAng_GCS = wrapTo180(SpinCalc('DCMtoEA123',legRR,.2,1));
legLAng_GCS = wrapTo180(SpinCalc('DCMtoEA123',legLR,.2,1));
shankRAng_GCS = wrapTo180(SpinCalc('DCMtoEA123',shankRR,.2,1)); 
shankLAng_GCS = wrapTo180(SpinCalc('DCMtoEA123',shankLR,.2,1));
footRAng_GCS = wrapTo180(SpinCalc('DCMtoEA123',footRR,.2,1));
footLAng_GCS = wrapTo180(SpinCalc('DCMtoEA123',footLR,.2,1));




%Calculate joint angles
% legRAng_pelvis = wrapTo180(SpinCalc('DCMtoEA123',((inv(pelvisR))*legRR),.2,1));
% legLAng_pelvis = wrapTo180(SpinCalc('DCMtoEA123',((inv(pelvisR))*legLR),.2,1));
% shankRAng_legR = wrapTo180(SpinCalc('DCMtoEA123',((inv(legRR))*shankRR),.2,1)); %precision tolerance violation 
% shankLAng_legL = wrapTo180(SpinCalc('DCMtoEA123',((inv(legLR))*shankLR),.2,1));
% footRAng_shankR = wrapTo180(SpinCalc('DCMtoEA123',((inv(shankRR))*footRR),.2,1)); %precision tolerance violation 
% footLAng_shankL = wrapTo180(SpinCalc('DCMtoEA123',((inv(shankLR))*footLR),.2,1));
% pelvisAng_GCS = wrapTo180(pelvisAng_GCS);

legRAng_pelvis = wrapTo180(SpinCalc('DCMtoEA123',((inv(pelvisR))*legRR),.2,1));
legLAng_pelvis = wrapTo180(SpinCalc('DCMtoEA123',((inv(pelvisR))*legLR),.2,1));
shankRAng_legR = wrapTo180(SpinCalc('DCMtoEA123',((inv(legRR))*shankRR),.2,1)); %precision tolerance violation 
shankLAng_legL = wrapTo180(SpinCalc('DCMtoEA123',((inv(legLR))*shankLR),.2,1));
footRAng_shankR = wrapTo180(SpinCalc('DCMtoEA123',((inv(shankRR))*footRR),.2,1)); %precision tolerance violation 
footLAng_shankL = wrapTo180(SpinCalc('DCMtoEA123',((inv(shankLR))*footLR),.2,1));
pelvisAng_GCS = wrapTo180(pelvisAng_GCS);

%Write segment and joint angles to data structure
static.GCSOrientations.Pelvis_GCS = pelvisAng_GCS;
static.GCSOrientations.LegLAng_GCS = legLAng_GCS;
static.GCSOrientations.LegRAng_GCS = legRAng_GCS;
static.GCSOrientations.ShankLAng_GCS = shankLAng_GCS;
static.GCSOrientations.ShankRAng_GCS = shankRAng_GCS;
static.GCSOrientations.footRAng_GCS = footRAng_GCS;
static.GCSOrientations.footLAng_GCS = footLAng_GCS;

static.km.Rankle = footRAng_shankR;
static.km.Lankle = footLAng_shankL;
static.km.Rhip = legRAng_pelvis;
static.km.Lhip = legLAng_pelvis;
static.km.Rknee = shankRAng_legR;
static.km.Lknee = shankLAng_legL;
static.km.PelvisAngle = pelvisAng_GCS;
static.km.pelvisCOMPos = static.markers.pelvisOrigin;

%% Process segment distal and proximal end positions and save to data structure (from v3d)
% 
% static.v3dOffsets.FAJCLPos_ShankL = FAJCLPos_ShankL;
% static.v3dOffsets.FAJCRPos_ShankR = FAJCRPos_ShankR;
% static.v3dOffsets.FAJCLPos_GCS = x*(mean(FAnkleL_landmark{1})'*100); 
% static.v3dOffsets.FAJCRPos_GCS = x*(mean(FAnkleR_landmark{1})'*100); clear FAJCLPos_ShankL FAJCRPos_ShankR FAnkleL_landmark FAnkleR_landmark
% 
% static.v3dOffsets.FHipLPos_LegL = FHipLPos_LegL;
% static.v3dOffsets.FHipRPos_LegR = FHipRPos_LegR;
% static.v3dOffsets.FHipLPos_GCS = x*(mean(FHipL_landmark{1})'*100); 
% static.v3dOffsets.FHipRPos_GCS = x*(mean(FHipR_landmark{1})'*100); clear FHipLPos_LegL FHipLPos_LegR FHipL_landmark FHipR_landmark
% 
% static.v3dOffsets.FKneeLPos_LegL = FKneeLPos_LegL;
% static.v3dOffsets.FKneeRPos_LegR = FKneeRPos_LegR;
% static.v3dOffsets.FKneeLPos_GCS = x*(mean(FKneeL_landmark{1})'*100); 
% static.v3dOffsets.FKneeRPos_GCS = x*(mean(FKneeR_landmark{1})'*100); clear FKneeLPos_LegL FKneeRPos_LegR F FKneeL_landmark FKneeR_landmark
% 
% static.v3dOffsets.LegLDistalPos_GCS = x*(mean(LegLDistalPos{1})'*100);
% static.v3dOffsets.LegLProxPos_GCS = x*(mean(LegLProxPos{1})'*100);
% static.v3dOffsets.LegRDistalPos_GCS = x*(mean(LegRDistalPos{1})'*100);
% static.v3dOffsets.LegRProxPos_GCS = x*(mean(LegRProxPos{1})'*100); 
% static.v3dOffsets.LegLProxPos_pelvis = x*(mean(LegLPos_pelvis{1})'*100);
% static.v3dOffsets.LegRProxPos_pelvis = x*(mean(LegRPos_pelvis{1})'*100);  
% static.v3dOffsets.ShankLDistalPos_GCS = x*(mean(ShankLDistalPos{1})'*100);
% static.v3dOffsets.ShankLProxPos_GCS = x*(mean(ShankLProxPos{1})'*100);
% static.v3dOffsets.ShankRDistalPos_GCS = x*(mean(ShankRDistalPos{1})'*100);
% static.v3dOffsets.ShankRProxPos_GCS = x*(mean(ShankRProxPos{1})'*100); 
% static.v3dOffsets.ShankLProxPos_LegL = x*(mean(ShankLPos_LegL{1})'*100);
% static.v3dOffsets.ShankRProxPos_LegR = x*(mean(ShankRPos_LegR{1})'*100);  
% static.v3dOffsets.PelvisDistPos_GCS = x*(mean(PelivsRDistalPos{1})'*100);
% static.v3dOffsets.PelvisProxPos_GCS = x*(mean(PelivsRProxPos{1})'*100);  

%% Estimate segment offsets
static.offsets.World2Pelvis = static.markers.pelvisOrigin;
static.offsets.Pelvis2HJCR = static.markers.HJCR - static.markers.pelvisOrigin;
static.offsets.Pelvis2HJCL = static.markers.HJCL - static.markers.pelvisOrigin;
static.offsets.HJCR2LegR =  static.markers.HJCR - (x*(mean(LegRProxPos{1})'*100));
static.offsets.HJCL2LegL =  static.markers.HJCL - (x*(mean(LegLProxPos{1})'*100));
% static.offsets.LegR2KJCR = static.markers.KJCR - (static.markers.HJCR + static.offsets.HJCR2LegR);
% static.offsets.LegL2KJCL = static.markers.KJCL - (static.markers.HJCL + static.offsets.HJCL2LegL);
static.offsets.HJCR2KJCR = static.markers.KJCR - static.markers.HJCR; %double check these
static.offsets.HJCL2KJCL = static.markers.KJCL - static.markers.HJCL;

static.offsets.KJCR2ShankR = static.markers.KJCR - (x*(mean(ShankRProxPos{1})'*100)) ;
static.offsets.KJCL2ShankL = static.markers.KJCL - (x*(mean(ShankLProxPos{1})'*100)) ;
static.offsets.ShankR2AJC_pros = static.markers.AJCR - (x*(mean(ShankRProxPos{1})'*100));
static.offsets.ShankL2AJCL = static.markers.AJCL - (x*(mean(ShankLProxPos{1})'*100));

%% Generate limb length estimates
%Estimate segment lengths
legLengthL = norm(static.markers.HJCL - static.markers.KJCL);
legLengthR = norm(static.markers.HJCR - static.markers.KJCR);

shankLengthL =  norm(static.markers.KJCL - static.markers.AJCL);
shankLengthR = norm(static.markers.KJCR(3) - (static.markers.LSKmid(3)));

pelvisXdist = norm(static.markers.ASISmidpt - static.markers.PSISmidpt);
pelvisYdist = norm(static.markers.RASI - static.markers.LASI);
pelvisZdist = (abs(static.markers.pelvisOrigin(3) - static.markers.HJCR(3))); 
ASISYdist = (abs(static.markers.RASI(2) - static.markers.LASI(2)));
PSISYdist = (abs(static.markers.RPSI(2) - static.markers.LPSI(2)));

footLdistX = abs(static.markers.LTCE(1) - static.markers.LHCE(1)); 
footLdistY = abs(static.markers.LTLA(2) - static.markers.LTME(2));
footLdistZ = static.markers.AJCL(3) - static.markers.ground(3);

% Save segment length estimates to data structure 
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

%% Filter force data
% Downsample forces and filter forces and COP
order = 2;
fs = 1000;
fc = 50;%50 hz low pass cut off
Wn_fp = 2*fc/fs; 
[b,a] = butter(order,Wn_fp); 
FP1filt = filtfilt(b,a,FP1{1}); %symmetric lowpass filter of normal force
FP2filt = filtfilt(b,a,FP2{1});
FP1COPfilt = filtfilt(b,a,FP1COP{1});
FP2COPfilt = filtfilt(b,a,FP2COP{1});
clear fc Wn_fp a b 

FP1ds = downsample(FP1filt,5);
FP2ds = downsample(FP2filt,5);
FP1COPds = downsample(FP1COPfilt,5);
FP2COPds = downsample(FP2COPfilt,5);
clear FP1filt FP2filt FP1COPfilt FP2COPfilt

%Rotate force data
FP1 = x*(mean(FP1ds)');
FP2 = x*(mean(FP2ds)');
FP1COP = x*(mean(FP1COPds)');
FP2COP = x*(mean(FP2COPds)');

%Save to data structure
static.forces.FP1 = FP1;
static.forces.FP2 = FP2;
static.forces.FP1COP = FP1COP;
static.forces.FP2COP = FP2COP;
static.framerate = FRAME_RATE;
static.trial = extractAfter(string(FILE_NAME{1}),"VSF\"); 
static.mass = mass{1};
static.height = height{1};
static.legLength = legLength;
clearvars -except static



