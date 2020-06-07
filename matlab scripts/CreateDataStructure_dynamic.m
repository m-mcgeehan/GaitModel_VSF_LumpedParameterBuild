%Code to:
% 1) read in subject kinematic data from visual 3d pipeline
% 2) align visual 3d coordinate system with simscape
% 3) filter kinematic and force data and downsample forces
% 4) create data structure to drive model (usable in "startup_GaitModel_VSF.m")

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

%% Create virtual markers and write them to the data structure
ASISmidpt = (markerSet.RASI+markerSet.LASI)/2;
PSISmidpt = (markerSet.RPSI+markerSet.LPSI)/2;
pelvisOrigin = PSISmidpt;
    ASISdist = norm(markerSet.RASI - markerSet.LASI);
HJCR = markerSet.FHipR_landmark;
HJCL = markerSet.FHipL_landmark;
KJCR = markerSet.FKneeR_landmark;
KJCL = markerSet.FKneeL_landmark;
FAJCL = markerSet.FAnkleL_landmark;
% HJCR = [markerSet.RASI(:,1) - (.36*ASISdist), markerSet.RASI(:,2) - (.14*ASISdist), markerSet.RASI(:,3) - (.22*ASISdist)];
% HJCL = [markerSet.LASI(:,1) - (.36*ASISdist), markerSet.LASI(:,2) + (.14*ASISdist), markerSet.LASI(:,3) - (.22*ASISdist)];
% KJCR = (markerSet.RLEP + markerSet.RMEP)/2;
% KJCL = (markerSet.LLEP + markerSet.LMEP)/2;
AJC_pros = [markerSet.RTCE(:,1) - markerSet.RHCE(:,1), markerSet.RHLA(:,2)+.4, markerSet.RHLA(:,3)+.62];
% AJCL = (markerSet.RHLA + markerSet.RHME)/2;

markerSet.ASISmidpt = ASISmidpt;
markerSet.PSISmidpt = PSISmidpt;
markerSet.pelvisOrigin = pelvisOrigin;
markerSet.HJCR = markerSet.FHipR_landmark;
markerSet.HJCL = markerSet.FHipL_landmark;
markerSet.KJCR = markerSet.FKneeR_landmark;
markerSet.KJCL = markerSet.FKneeL_landmark;
markerSet.AJCL = FAJCL;
markerSet.AJC_pros = AJC_pros;
markerSet.FAJC_pros = markerSet.FAnkleR_landmark;

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

%% Estimate joint angle data and create data structure
numFrames = length(time);
%pelvis
for i = 1:numFrames
    pelvis_j = (markerSet.LPSI(i,:) - markerSet.pelvisOrigin(i,:))/norm(markerSet.LPSI(i,:)-markerSet.pelvisOrigin(i,:));
    pelvis_k = cross((markerSet.ASISmidpt(i,:) - markerSet.PSISmidpt(i,:))/(norm(markerSet.ASISmidpt(i,:)-markerSet.PSISmidpt(i,:))),pelvis_j);
    pelvis_i = cross(pelvis_j,pelvis_k);
    pelvisR{i} = [pelvis_i', pelvis_j', pelvis_k'];
    if det(pelvisR{i}) < .9
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
    if det(legRR{i}) < .9
    warning('determinant of legRR is below threshold. Determinant = ');
    display(det(legRR{i}));
    else
    end
end

%LegL
for i = 1:numFrames
    legL_k = (markerSet.HJCL(i,:) - markerSet.KJCL(i,:))/(norm(markerSet.HJCL(i,:) - markerSet.KJCL(i,:)));
    legL_i = cross(((markerSet.LLEP(i,:) - markerSet.LMEP(i,:))/norm(markerSet.LLEP(i,:) - markerSet.LMEP(i,:))),legL_k);
    legL_j = cross(legL_k,legL_i);
    legLR{i} = [legL_i', legL_j', legL_k'];
    if det(legLR{i}) < .9
    warning('determinant of legLR is below threshold. Determinant = ');
    display(det(legLR{i}));
    else
    end
end

%ShankR
for i = 1:numFrames
    shankR_k = (markerSet.KJCR(i,:) - markerSet.FAJC_pros(i,:))/(norm(markerSet.KJCR(i,:) - markerSet.FAJC_pros(i,:)));
    shankR_i = cross(((markerSet.RMEP(i,:) - markerSet.RLEP(i,:))/norm(markerSet.RMEP(i,:) - markerSet.RLEP(i,:))),shankR_k);
    shankR_j = cross(shankR_k,shankR_i);
    shankRR{i} = [shankR_i', shankR_j', shankR_k'];
    if det(shankRR{i}) < .9
    warning('determinant of shankRR is below threshold. Determinant = ');
    display(det(shankRR{i}));
    else
    end
end

%ShankL
for i = 1:numFrames
   shankL_k = (markerSet.KJCL(i,:) - markerSet.AJCL(i,:))/(norm(markerSet.KJCL(i,:) - markerSet.AJCL(i,:)));
   shankL_i = cross(((markerSet.LLEP(i,:) - markerSet.LMEP(i,:))/norm(markerSet.LLEP(i,:) - markerSet.LMEP(i,:))),shankL_k);
   shankL_j = cross(shankL_k,shankL_i);
   shankLR{i} = [shankL_i', shankL_j', shankL_k'];
   if det(shankLR{i}) < .9
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
    if det(footLR{i}) < .9
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
    LAnkleAngle = LAnkleAngle;
    RKneeAngle = RKneeAngle;
    LKneeAngle = LKneeAngle;
    RHipAngle = RHipAngle;
    LHipAngle = LHipAngle;
end

    
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
fc = 50;%50 hz low pass cut off
Wn_fp = 2*fc/fs; 
[b,a] = butter(order,Wn_fp); 
FP1filt = filtfilt(b,a,FP1{1}); %symmetric lowpass filter of forces
FP2filt = filtfilt(b,a,FP2{1});

clear fc Wn_fp a b

FP1ds = downsample(FP1filt,5);
FP2ds = downsample(FP2filt,5);

%% set up structure
subj1.km.Lankle = LAnkleAngle;
subj1.km.Rankle = RAnkleAngle;
subj1.km.Rhip = RHipAngle;
subj1.km.Lhip = LHipAngle;
subj1.km.Rknee = RKneeAngle;
subj1.km.Lknee = LKneeAngle;
subj1.km.PelvisAngle = pelvisAng_GCS;
subj1.km.PelvisCOMPos = pelvisOrigin;
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
