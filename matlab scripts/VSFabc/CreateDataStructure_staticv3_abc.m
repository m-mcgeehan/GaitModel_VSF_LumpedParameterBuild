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

%% Load mat file
[fName,fPath] = uigetfile('.mat');
load(fullfile(fPath,fName));
%% Rotate markers and setup marker structure
x = [0 -1 0;...  %rotation matrix to align coordinate systems (simscape wrt v3d: x = y, y = -x, z = z)
    1 0 0;...
    0 0 1];

% markerNames = {'FHipR_landmark' 'FHipL_landmark' 'FKneeL_landmark' 'FKneeR_landmark' 'FAnkleL_landmark' 'FAnkleR_landmark' 'LASI' 'RASI' 'LPSI'...
%     'RPSI' 'LGTR' 'RGTR' 'LTH2' 'LTH3' 'LTH4' 'RTH1' 'RTH2' 'RTH3' 'RTH4' 'LMEP' 'LLEP' 'RMEP' 'RLEP' 'LSK1' 'LSK2' 'LSK3' 'LSK4'...
%      'LHCE' 'RHCE' 'LHLA' 'RHLA' 'LHME' 'RHME' 'LTCE' 'RTCE' 'LTLA' 'RTLA' 'LTME' 'RTME'};

markerNames = {'RTME' 'RTLA' 'RTH4' 'RTH3' 'RTH2' 'RTH1' 'RTCE' 'RPSI' 'RMEP' 'RLEP' 'RHME' 'RHLA' 'RHCE' 'RGTR' 'RASI' 'LTME'...
    'LTLA' 'LTH4' 'LTH3' 'LTH2' 'LTH1' 'LTCE' 'LSK4' 'LSK3' 'LSK2' 'LSK1' 'LPSI' 'LMEP' 'LLEP' 'LHME' 'LHLA' 'LHCE' 'LGTR' 'LASI'...
    'RKJC' 'RHJC_reg' 'RHJC_fxnl' 'RAJC_fxnl' 'LKJC' 'LHJC_reg' 'LHJC_fxnl' 'LAJC_fxnl' 'VSF_Rankle' 'CODA_origin'};

for i = 1:length(markerNames) %Scan through marker list, rotate, convert from m to cm, write to data structure
    marker = markerNames{i};
    eval(['markerData =',marker,';']);
    markerDataRot = 100*(x*mean(markerData{1}(:,1:3))');
    static.markers.(marker) = markerDataRot;
end

static.markers.ground = [0 0 0];

clear 'RTME' 'RTLA' 'RTH4' 'RTH3' 'RTH2' 'RTH1' 'RTCE' 'RPSI' 'RMEP' 'RLEP' 'RHME' 'RHLA' 'RHCE' 'RGTR' 'RASI' 'LTME'...
    'LTLA' 'LTH4' 'LTH3' 'LTH2' 'LTH1' 'LTCE' 'LSK4' 'LSK3' 'LSK2' 'LSK1' 'LMEP' 'LLEP' 'LHME' 'LHLA' 'LHCE' 'LGTR' 'LASI'...
    'RKJC' 'RHJC_reg' 'RHJC_fxnl' 'RAJC_fxnl' 'LKJC' 'LHJC_reg' 'LHJC_fxnl' 'LAJC_fxnl' 'VSF_Rankle' 'CODA_origin'
%%
%format virtual markers
ASISmidpt = (static.markers.RASI+static.markers.LASI)/2;
PSISmidpt = (static.markers.RPSI+static.markers.LPSI)/2;
pelvisOrigin = ASISmidpt;
    ASISdist = norm(static.markers.RASI - static.markers.LASI);
% HJCR = (x*(mean(FHipR_landmark{1}(:,1:3))')*100);
% HJCL = (x*(mean(FHipL_landmark{1}(:,1:3))')*100);
legLength = 10*(norm(ASISmidpt - static.markers.LMEP) + norm(static.markers.LMEP - static.markers.LHME));
HJCR_reg = pelvisOrigin + (.1*[(11-(.063*legLength)); -1*(8+(.086*legLength)); (9-(.078*legLength))]);
HJCL_reg = pelvisOrigin + (.1*[(11-(.063*legLength)); (8+(.086*legLength)); (9-(.078*legLength))]);
KJCR_midpt = (static.markers.RMEP + static.markers.RLEP)/2;
KJCL_midpt = (static.markers.LMEP + static.markers.LLEP)/2;
% KJCR = (x*(mean(FKneeR_landmark{1}(:,1:3))')*100);
% KJCL = (x*(mean(FKneeL_landmark{1}(:,1:3))')*100);
AJC_pros = [static.markers.RTCE(1) - static.markers.RHCE(1), static.markers.RHLA(2)+.4, static.markers.RHLA(3)+.62]';
AJCL = (static.markers.LHLA + static.markers.LHME)/2;
pelvisAng_GCS_v3d = pelvisAng_GCS;

static.markers.ASISmidpt = ASISmidpt;
static.markers.PSISmidpt = PSISmidpt;
static.markers.pelvisOrigin = pelvisOrigin;
static.markers.HJCR = HJCR_reg;
static.markers.HJCL = HJCL_reg;
static.markers.KJCR = KJCR_midpt;
static.markers.KJCL = KJCL_midpt;
static.markers.AJC_pros = AJC_pros;
static.markers.AJCR = (x*(mean(FAnkleR_landmark{1}(:,1:3))')*100);
% static.markers.AJCL = (x*(mean(FAnkleL_landmark{1}(:,1:3))')*100);
static.markers.AJCL = AJCL;
clear ASISmidpt PSISmidpt pelvisOrigin HJCR HJCL KJCR KJCL AJC_pros AJCR AJCL 
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
% legR_k = (static.markers.HJCR - static.markers.KJCR)/(norm(static.markers.HJCR - static.markers.KJCR));
% legR_i = cross(((static.markers.RMEP - static.markers.RLEP)/norm(static.markers.RMEP - static.markers.RLEP)),legR_k);
% legR_j = cross(legR_k,legR_i);
% legRR = [legR_i, legR_j, legR_k]';
% 
% if det(legRR) < .95
% warning('determinant of legRR is below threshold. Determinant = ');
% display(det(legRR));
% else
% end

legR_k = (static.markers.RHJC_fxnl - static.markers.RKJC)/(norm(static.markers.RHJC_fxnl - static.markers.RKJC));
legR_i = cross(((static.markers.RMEP - static.markers.RLEP)/norm(static.markers.RMEP - static.markers.RLEP)),legR_k);
legR_j = cross(legR_k,legR_i);
legRR = [legR_i, legR_j, legR_k]';

if det(legRR) < .95
warning('determinant of legRR is below threshold. Determinant = ');
display(det(legRR));
else
end

% legL_k = (static.markers.HJCL - static.markers.KJCL)/(norm(static.markers.HJCL - static.markers.KJCL));
% legL_i = cross(((static.markers.LLEP - static.markers.LMEP)/norm(static.markers.LLEP - static.markers.LMEP)),legL_k); 
% legL_j = cross(legL_k,legL_i);
% legLR = [legL_i, legL_j, legL_k]';
% legLRdet = det(legLR);
% 
% if det(legLR) < .95
% warning('determinant of legLR is below threshold. Determinant = ');
% display(det(legLR));
% else
% end
% legL_k = (static.markers.LHJC_fxnl - static.markers.LKJC)/(norm(static.markers.LHJC_fxnl - static.markers.LKJC));
% legL_i = cross(((static.markers.LLEP - static.markers.LMEP)/norm(static.markers.LLEP - static.markers.LMEP)),legL_k); 
% legL_j = cross(legL_k,legL_i);

legL_k = (static.markers.LHJC_reg - static.markers.LKJC)/(norm(static.markers.LHJC_reg - static.markers.LKJC));
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
% shankR_k = (static.markers.KJCR - static.markers.AJCR)/(norm(static.markers.KJCR - static.markers.AJCR));
% shankR_i = cross(((static.markers.RMEP - static.markers.RLEP)/norm(static.markers.RMEP - static.markers.RLEP)),shankR_k);
% shankR_j = cross(shankR_k,shankR_i);
% shankRR = [shankR_i, shankR_j, shankR_k]';
% 
% if det(shankRR) < .95
% warning('determinant of shankRR is below threshold. Determinant = ');
% display(det(shankRR));
% else
% end

shankR_k = (static.markers.RKJC - static.markers.VSF_Rankle)/(norm(static.markers.RKJC - static.markers.VSF_Rankle));
shankR_i = cross(((static.markers.RMEP - static.markers.RLEP)/norm(static.markers.RMEP - static.markers.RLEP)),shankR_k);
shankR_j = cross(shankR_k,shankR_i);
shankRR = [shankR_i, shankR_j, shankR_k]';

if det(shankRR) < .95
warning('determinant of shankRR is below threshold. Determinant = ');
display(det(shankRR));
else
end
% 
% shankL_k = (static.markers.KJCL - static.markers.AJCL)/(norm(static.markers.KJCL - static.markers.AJCL));
% shankL_i = cross(((static.markers.LLEP - static.markers.LMEP)/norm(static.markers.LLEP - static.markers.LMEP)),shankL_k);
% shankL_j = cross(shankL_k,shankL_i);
% shankLR = [shankL_i, shankL_j, shankL_k];
% 
% if det(shankLR) < .95
% warning('determinant of shankLR is below threshold. Determinant = ');
% display(det(shankLR));
% else
% end


% shankL_k = (static.markers.LKJC - static.markers.LAJC_fxnl)/(norm(static.markers.LKJC - static.markers.LAJC_fxnl));
% shankL_i = cross(((static.markers.LLEP - static.markers.LMEP)/norm(static.markers.LLEP - static.markers.LMEP)),shankL_k);
% shankL_j = cross(shankL_k,shankL_i);
% shankLR = [shankL_i, shankL_j, shankL_k];

shankL_k = (static.markers.LKJC - static.markers.AJCL)/(norm(static.markers.LKJC - static.markers.AJCL));
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

%% V3d angles


jointAngleNames = {'LAnkleAngle' 'LHipAngle' 'LKneeAngle' 'RAnkleAngle' 'RHipAngle' 'RKneeAngle' 'pelvisAng_GCS_v3d'};

for i = 1:length(jointAngleNames) %Scan through marker list, rotate, convert from m to cm, write to data structure
    joint = jointAngleNames{i};
    eval(['jointAngle =',joint,';']);
    jointAngleRot = (x*mean(jointAngle{1}(:,1:3))');
    static.km.v3d.(joint) = jointAngleRot;
end

%% Process segment distal and proximal end positions and save to data structure (from v3d)
offsetNames = {'FAJCLPos_ShankL' 'FAJCRPos_ShankR' 'FHipLPos_LegL' 'FHipRPos_LegR' 'FKneeLPos_LegL' 'FKneeRPos_LegR' 'FKneeRPos_ShankR'...
    'FKneeLPos_ShankL' 'FootLPos_ShankL' 'FootRPos_ShankR' 'LegLPos_pelvis' 'LegRPos_pelvis' 'ShankRPos_LegR' 'ShankLPos_LegL'...
    'FootLPos_ShankL' 'FootRPos_ShankR'};

for i = 1:length(offsetNames) %Scan through offset list, rotate, convert from m to cm, write to data structure
    offset = offsetNames{i};
    eval(['offsetData =',offset,';']);
    offsetDataRot = 100*(x*mean(offsetData{1}(:,1:3))');
    static.offsets.(offset) = offsetDataRot;
end

clear 'FAJCLPos_ShankL' 'FAJCRPos_ShankR''FHipLPos_LegL' 'FHipRPos_LegR' 'FKneeLPos_LegL' 'FKneeRPos_LegR' 'FKneeRPos_ShankR'...
    'FKneeLPos_ShankL' 'FootLPos_ShankL' 'FootRPos_ShankR' 'LegLPos_pelvis' 'LegRPos_pelvis' 'ShankRPos_LegR' 'ShankLPos_LegL'...
    'FootLPos_ShankL' 'FootRPos_ShankR'

oTrans = static.markers.pelvisOrigin - static.markers.PSISmidpt;
static.offsets.LegLPos_pelvis1 = static.offsets.LegLPos_pelvis - oTrans;
static.offsets.LegRPos_pelvis1 = static.offsets.LegRPos_pelvis - oTrans;

%% Estimate segment offsets from marker data


static.offsets.World2Pelvis = static.markers.pelvisOrigin;
static.offsets.Pelvis2HJCR = static.offsets.LegRPos_pelvis;
static.offsets.Pelvis2HJCL = static.offsets.LegLPos_pelvis;
static.offsets.HJCR2LegR =  static.offsets.FHipRPos_LegR;
static.offsets.HJCL2LegL =  static.offsets.FHipLPos_LegL;

static.offsets.HJCR2KJCR = static.offsets.FKneeRPos_LegR; %double check these
static.offsets.HJCL2KJCL = static.offsets.FKneeLPos_LegL;

static.offsets.KJCR2ShankR = static.offsets.FKneeRPos_ShankR;
static.offsets.KJCL2ShankL = static.offsets.FKneeLPos_ShankL;

static.offsets.ShankR2AJC_pros = static.offsets.FAJCRPos_ShankR;
static.offsets.ShankL2AJCL = static.offsets.FAJCLPos_ShankL;

% static.offsets.World2Pelvis = static.markers.pelvisOrigin;
% static.offsets.Pelvis2HJCR = static.markers.HJCR - static.markers.pelvisOrigin;
% static.offsets.Pelvis2HJCL = static.markers.HJCL - static.markers.pelvisOrigin;
% static.offsets.HJCR2LegR =  static.markers.HJCR - (x*(mean(LegRProxPos{1})'*100));
% static.offsets.HJCL2LegL =  static.markers.HJCL - (x*(mean(LegLProxPos{1})'*100));
% % static.offsets.LegR2KJCR = static.markers.KJCR - (static.markers.HJCR + static.offsets.HJCR2LegR);
% % static.offsets.LegL2KJCL = static.markers.KJCL - (static.markers.HJCL + static.offsets.HJCL2LegL);
% static.offsets.HJCR2KJCR = static.markers.KJCR - static.markers.HJCR; %double check these
% static.offsets.HJCL2KJCL = static.markers.KJCL - static.markers.HJCL;
% 
% static.offsets.KJCR2ShankR = static.markers.KJCR - (x*(mean(ShankRProxPos{1})'*100)) ;
% static.offsets.KJCL2ShankL = static.markers.KJCL - (x*(mean(ShankLProxPos{1})'*100)) ;
% static.offsets.ShankR2AJC_pros = static.markers.AJCR - (x*(mean(ShankRProxPos{1})'*100));
% static.offsets.ShankL2AJCL = static.markers.AJCL - (x*(mean(ShankLProxPos{1})'*100));



%% Generate limb length estimates
%Estimate segment lengths
% legLengthL = norm(static.markers.HJCL - static.markers.KJCL);
% legLengthR = norm(static.markers.HJCR - static.markers.KJCR);
% 
% shankLengthL =  norm(static.markers.KJCL - static.markers.AJCL);
% shankLengthR = norm(static.markers.KJCR(3) - (static.markers.LSKmid(3)));

legLengthL = LTH_Length{1};
legLengthR = RTH_Length{1};

shankLengthL = LSK_Length{1};
shankLengthR = RSK_Length{1}; clear LTH_Length RTH_Length  LSK_Length  RSK_Length


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
static.trial = fName;
static.mass = mass{1};
static.height = height{1};
static.legLength = legLength;
clearvars -except static

%Name and save variable

varName = convertStringsToChars(extractBefore(static.trial,"."));

eval([(varName) '=' 'static;']); 

% eval(['(convertStringsToChars(varName))' ' = static;']);

% [file,path,indx] = uiputfile(strcat(varName,'.m'));

uisave((varName),(varName));
clear static file varName path indx

