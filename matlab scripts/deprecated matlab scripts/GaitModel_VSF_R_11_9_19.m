%% joint rotations
Rleft = [-1 0 0;...
        0 -1 0;...
        0 0 -1];
R = [1 0 0;...
    0 -1 0;...
    0 0 1];
x = [0 1 0;...
    -1 0 0;...
    0 0 1];

%Pelvis
ASISmid = (subj1.static.markers.RASI+subj1.static.markers.LASI)/2;
PSISmid = (subj1.static.markers.RPSI+subj1.static.markers.LPSI)/2;
pelvisOrigin = (ASISmid+PSISmid)/2;
PelvisAng_GCS = R*subj1.static.GCSOrientations.Pelvis_GCS;%Rotation %%NEED TO ZERO IN THE CAD MODEL
PelvisPos_GCS = pelvisOrigin; %Translation

%Legs
FHJCRPos_pelvis = subj1.static.offsets.FHipRPos_GCS - PelvisPos_GCS ; %Translation 
FHJCLPos_pelvis = subj1.static.offsets.FHipLPos_GCS - PelvisPos_GCS ; %Translation 
LegRAng_pelvis = PelvisAng_GCS - R*subj1.static.GCSOrientations.LegRAng_GCS; %Rotation
LegLAng_pelvis = PelvisAng_GCS - Rleft*subj1.static.GCSOrientations.LegLAng_GCS;%Rotation
LegRPos_HJCR = subj1.static.offsets.LegRProxPos_GCS - subj1.static.offsets.FHipRPos_GCS;
LegLPos_HJCR = subj1.static.offsets.LegLProxPos_GCS - subj1.static.offsets.FHipLPos_GCS;
FKJCRPos_LegR = ((subj1.static.markers.RMEP + subj1.static.markers.RLEP)/2) - subj1.static.offsets.LegRProxPos_GCS; %x is opposite
FKJCLPos_LegL = (subj1.static.markers.LMEP + subj1.static.markers.LLEP)/2 - subj1.static.offsets.LegLProxPos_GCS;
% FKJCRPos_LegR = subj1.static.offsets.FKneeLPos_LegL;
% FKJCLPos_LegL = subj1.static.offsets.FKneeRPos_GCS - subj1.static.offsets.LegRProxPos_GCS; %%Lateral displacement

%Shanks
ShankRAng_LegR = R*subj1.static.GCSOrientations.LegRAng_GCS - R*subj1.static.GCSOrientations.ShankRAng_GCS;
ShankLAng_LegL = Rleft*subj1.static.GCSOrientations.LegLAng_GCS - Rleft*subj1.static.GCSOrientations.ShankLAng_GCS;
ShankRPos_KJCR = subj1.static.offsets.ShankRProxPos_GCS - ((subj1.static.markers.RMEP + subj1.static.markers.RLEP)/2);
ShankLPos_KJCL = subj1.static.offsets.ShankLProxPos_GCS - (subj1.static.markers.LMEP + subj1.static.markers.LLEP)/2;
FAJCLPos_ShankL = subj1.static.offsets.FAJCLPos_GCS - subj1.static.offsets.ShankLProxPos_GCS;
socket2shank_offset = [1.5 1.5 subj1.static.segLength.shankLengthR+3.5]; %cm 
pylon_length = (subj1.static.segLength.shankLengthL - subj1.static.segLength.shankLengthR) - 3.5 - (subj1.static.segLength.footLZ - 9.3);
%Feet
FootLAng_ShankL = Rleft*subj1.static.km.Lankle;

%Legs
RLegWrtPelvis = R*(subj1.static.km.Rhip);%rotation
LLegWrtPelvis = Rleft*(subj1.static.km.Lhip);%rotation
RLeg2PelvisOffset = (subj1.static.offsets.PelvisProxPos_GCS - subj1.static.offsets.LegRProxPos_GCS); %translation
LLeg2PelvisOffset = (subj1.static.offsets.PelvisProxPos_GCS - subj1.static.offsets.LegLProxPos_GCS); %translation

LegR2KJCOffset = R*subj1.static.offsets.FKneeRPos_LegR;%translation
LegL2KJCOffset=  Rleft*subj1.static.offsets.FKneeLPos_LegL;%translation

ShankRWrtLegR = R*subj1.static.km.Rknee;%rotation
ShankLWrtLegL = Rleft*subj1.static.km.Lknee;%rotation
ShankR2KJCOffset = (subj1.static.offsets.FKneeRPos_GCS - subj1.static.offsets.ShankRProxPos_GCS);%translation
ShankL2KJCOffset = (subj1.static.offsets.FKneeLPos_GCS - subj1.static.offsets.ShankLProxPos_GCS);%translation

FootLWrtShankL = R*subj1.static.km.Rknee;%rotation
ShankL2AJCOffset = subj1.static.offsets.ShankLProxPos_GCS - subj1.static.offsets.FAJCLPos_GCS;%translation



%% TEMPS
%RIGHT SIDE
x = [0 0 0];
RLegWrtPelvis = x;%upper leg to hip
RLeg2KJCOffset =x;% knee to upper leg
RShank2Knee_offset = x;% lower leg to knee
RShankWrtLeg = x;

%Left SIDE
x = [0 0 0];
LLegWrtPelvis = x;%upper leg to hip
LLeg2KJCOffset =x;% knee to upper leg
LShank2Knee_offset = x;% lower leg to knee
LShankWrtLeg = x;
LShank2Ankle_offset = x; %ankle to lower leg

% subj1.offsets.RLegwrtPelvis;%translation
% LLeg2PelvisOffset = subj1.offsets.LLegwrtPelvis;%translation
   z = [0 -1 0;...
        1 0 0;...
        0 0 -1];
    y = [-1 0 0;...
        0 -1 0;...
        0 0 1];
RLeg2KJCOffset = subj1.static.offsets.FKneeRPos_LegR;%translation ***THIS NEEDS TO BE THE CENTER MARKER from v3d
LLeg2KJCOffset = subj1.static.offsets.FKneeLPos_LegL;%translation ***THIS NEEDS TO BE THE CENTER MARKER from v3d
%Shanks


RShankWrtLeg = R*(mean(subj1.static.km.Rknee)');
LShankWrtLeg = Rleft*(mean(subj1.static.km.Lknee)');


RHip2Torso_offset = subj1.segOffset.RLegwrtPelvis; %cm
LHip2Torso_offset = subj1.segOffset.LLegwrtPelvis;
    y = [1 0 0;...
        0 1 0;...
        0 0 -1];
    z = [1 0 0;...
        0 1 0;...
        0 0 -1];
RKnee2Leg_offset = y*subj1.segOffset.RShankwrtLeg; clear y
LKnee2Leg_offset = z*subj1.segOffset.LShankwrtLeg; clear z
LShank2Knee_offset = osf*[0 0 2];
RShank2Knee_offset = osf*[0 0 2];
LShank2Ankle_offset = osf*[2 0 subj1.seg.shankLengthL+.75];

init_height = subj1.static.offsets.PelvisProxPos_GCS(3);

marker_color = [0.0 0.4 1.0];
marker_radius = 0.5; %cm