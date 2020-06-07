% load('s1static');
static = s1static;
dynamic = subj1;

%Average static trial angles
PelvisStatic = [mean(static.km.PelvisAngle(:,1)),mean(static.km.PelvisAngle(:,2)),...
    mean(static.km.PelvisAngle(:,3))];
RhipStatic = [mean(static.km.Rhip(:,1)),mean(static.km.Rhip(:,2)),...
    mean(static.km.Rhip(:,3))];
LhipStatic = [mean(static.km.Lhip(:,1)),mean(static.km.Lhip(:,2)),...
    mean(static.km.Lhip(:,3))];
RkneeStatic = [mean(static.km.Rknee(:,1)),mean(static.km.Rknee(:,2)),...
    mean(static.km.Rknee(:,3))];
LkneeStatic = [mean(static.km.Lknee(:,1)),mean(static.km.Lknee(:,2)),...
    mean(static.km.Lknee(:,3))];
RankleStatic = [mean(static.km.Rankle(:,1)),mean(static.km.Rankle(:,2)),...
    mean(static.km.Rankle(:,3))];
LankleStatic = [mean(static.km.Lankle(:,1)),mean(static.km.Lankle(:,2)),...
    mean(static.km.Lankle(:,3))];

%Normalize dynamic trials to static values
PelvisAngNorm = dynamic.km.PelvisAngle - PelvisStatic;
RhipNorm = dynamic.km.Rhip - RhipStatic;
LhipNorm = dynamic.km.Lhip - LhipStatic;
RkneeNorm = dynamic.km.Rknee - RkneeStatic;
LkneeNorm = dynamic.km.Lknee - LkneeStatic;
% RankleNorm = dynamic.km.Rankle - RankleStatic;
LankleNorm = dynamic.km.Lankle - LankleStatic;


% PelvisAngNorm = [dynamic.km.PelvisAngle(:,1),dynamic.km.PelvisAngle(:,2) - PelvisStatic(:,2),dynamic.km.PelvisAngle(:,3)];
%RankleNorm = [linspace(0,0,length(subj1norm.km.Rhip(:,1)));linspace(0,0,length(subj1norm.km.Rhip(:,1)));linspace(0,0,length(subj1norm.km.Rhip(:,1)))];
%LankleNorm = [linspace(0,0.001,length(dynamic.km.PelvisAngle(:,1)));linspace(0,0.001,length(dynamic.km.PelvisAngle(:,1)));linspace(0,0.001,length(dynamic.km.PelvisAngle(:,1)))];
RankleNorm = [linspace(0,0,length(dynamic.km.PelvisAngle(:,1)));linspace(0,0,length(dynamic.km.PelvisAngle(:,1)));linspace(0,0,length(dynamic.km.PelvisAngle(:,1)))];


subj1norm.km.Rankle = RankleNorm';
subj1norm.km.Lankle = LankleNorm;
subj1norm.km.Rhip = RhipNorm;
subj1norm.km.Lhip = LhipNorm;
subj1norm.km.Rknee = dynamic.km.Rknee;
subj1norm.km.Lknee = dynamic.km.Lknee;
subj1norm.km.PelvisAngle = PelvisAngNorm;
subj1norm.km.PelvisCOMPos = dynamic.km.PelvisCOMPos;
subj1norm.forces.FP1 = dynamic.forces.FP1;
subj1norm.forces.FP2 = dynamic.forces.FP2;
subj1norm.framerate = static.framerate;
subj1norm.seg = dynamic.seg;
subj1norm.trial = subj1.trial;
% 
% subj1.km.Rankle = dynamic.km.Rankle;
% subj1.km.Lankle = dynamic.km.Lankle;
% subj1.km.Rhip = dynamic.km.Rhip;
% subj1.km.Lhip = dynamic.km.Lhip;
% subj1.km.Rknee = dynamic.km.Rknee;
% subj1.km.Lknee = dynamic.km.Lknee;
% subj1.km.PelvisAngle = PelvisAngNorm;
% subj1.km.PelvisCOMPos = dynamic.km.PelvisCOMPos;
% subj1.forces.FP1 = dynamic.forces.FP1;
% subj1.forces.FP2 = dynamic.forces.FP2;
% subj1.framerate = dynamic.framerate;
% 

% RankleNorm = RankleNorm';
% 
% subj1norm.km.Rankle = RankleNorm';
% subj1norm.km.Lankle = LankleNorm;
% subj1norm.km.Rhip = dynamic.km.Rhip;
% subj1norm.km.Lhip = dynamic.km.Lhip;
% subj1norm.km.Rknee = dynamic.km.Rknee;
% subj1norm.km.Lknee = dynamic.km.Lknee;
% subj1norm.km.PelvisAngle = PelvisAngNorm;
% subj1norm.km.PelvisCOMPos = dynamic.km.PelvisCOMPos;
% subj1norm.forces.FP1 = dynamic.forces.FP1;
% subj1norm.forces.FP2 = dynamic.forces.FP2;
% subj1norm.framerate = dynamic.framerate;

subj1 = subj1norm;
%% Figures...might be helpful
figure;
subplot(1,3,1);plot(PelvisAngNorm(:,1)); hold on; plot(dynamic.km.PelvisAngle(:,1));
subplot(1,3,2);plot(PelvisAngNorm(:,2)); hold on; plot(dynamic.km.PelvisAngle(:,2));
subplot(1,3,3);plot(PelvisAngNorm(:,3)); hold on; plot(dynamic.km.PelvisAngle(:,3));


figure;  
subplot(3,4,1);plot((subj1norm.km.Rhip(:,2)),'r');hold on;plot(dynamic.km.Rhip(:,2),'k');
title('Right Hip angleY')
subplot(3,4,2);plot((subj1norm.km.Lhip(:,2)),'r');hold on;plot(dynamic.km.Lhip(:,2),'k');
title('Left Hip angleY')
subplot(3,4,3);plot((subj1norm.km.Rknee(:,2)),'r');hold on;plot(dynamic.km.Rknee(:,2),'k');
title('Right Knee angleY')
subplot(3,4,4);plot((subj1norm.km.Lknee(:,2)),'r');hold on;plot(dynamic.km.Lknee(:,2),'k');
title('Left Knee angleY')
subplot(3,4,5);plot((subj1norm.km.Rankle(:,2)),'r');hold on;plot(dynamic.km.Rankle(:,2),'k');
title('Right Ankle angleY')
subplot(3,4,6);plot((subj1norm.km.Lankle(:,2)),'r');hold on;plot(dynamic.km.Lankle(:,2),'k');
title('Left Ankle angleY')
subplot(3,4,7);plot((subj1norm.km.PelvisAngle(:,2)),'r');hold on;plot(dynamic.km.PelvisAngle(:,2),'k');
title('Pelvis angleY')
subplot(3,4,8);plot((subj1norm.km.PelvisAngle(:,1)),'r');hold on;plot(dynamic.km.PelvisAngle(:,1),'k');
title('Pelvis angleX')
subplot(3,4,9);plot((subj1norm.km.PelvisAngle(:,3)),'r');hold on;plot(dynamic.km.PelvisAngle(:,3),'k');
title('Pelvis angleZ')
subplot(3,4,10);plot((dynamic.km.PelvisCOMPos(:,1)),'r');hold on;plot(dynamic.km.PelvisCOMPos(:,1),'k');
title('Pelvis PositionX')
sgtitle('norm (red) vs unnorm (black))');




figure;  
subplot(3,2,1);plot(time,static.km.Rhip(:,2),'k');
title('Right Hip angle')
subplot(3,2,2);plot(time,static.km.Lhip(:,2),'k');
title('Left Hip angle')
subplot(3,2,3);plot(time,static.km.Rknee(:,2),'k');
title('Right Knee angle')
subplot(3,2,4);plot(time,static.km.Lknee(:,2),'k');
title('Left Knee angle')
subplot(3,2,5);plot(time,static.km.Rankle(:,2),'k');
title('Right Ankle angle')
subplot(3,2,6);plot(time,static.km.Lankle(:,2),'k');
title('Left Ankle angle')
sgtitle('static trial sagittal plane kinematics)');

