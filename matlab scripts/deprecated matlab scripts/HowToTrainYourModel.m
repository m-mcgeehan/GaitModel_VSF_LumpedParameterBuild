%% Save sim results
simResults.pelvis = measPelvisPos;
simResults.R = measR;
simResults.L =measL;
simResults.fp = linspace(66,66,length(subj1.km.Rankle));

%% Create regression variables

%Concatenate data records for trials 1 and 2
s1conc.pelvis.px = [simResultsT1.pelvis.px.data; simResultsT2.pelvis.px.data];
s1conc.pelvis.vx = [simResultsT1.pelvis.vx.data; simResultsT2.pelvis.vx.data];
s1conc.pelvis.py = [simResultsT1.pelvis.py.data; simResultsT2.pelvis.py.data];
s1conc.pelvis.vy = [simResultsT1.pelvis.vy.data; simResultsT2.pelvis.vy.data];
s1conc.pelvis.pz = [simResultsT1.pelvis.pz.data; simResultsT2.pelvis.pz.data];
s1conc.pelvis.vz = [simResultsT1.pelvis.vz.data; simResultsT2.pelvis.vz.data];
s1conc.pelvis.qx = [simResultsT1.pelvis.qx.data; simResultsT2.pelvis.qx.data];
s1conc.pelvis.qy = [simResultsT1.pelvis.qy.data; simResultsT2.pelvis.qy.data];
s1conc.pelvis.qz = [simResultsT1.pelvis.qz.data; simResultsT2.pelvis.qz.data];

s1conc.RhipAngle = [simResultsT1.R.hip_angle.data; simResultsT2.R.hip_angle.data];
s1conc.RhipSpeed = [simResultsT1.R.hip_speed.data; simResultsT2.R.hip_speed.data];
s1conc.RkneeAngle = [simResultsT1.R.knee_angle.data; simResultsT2.R.knee_angle.data];
s1conc.RkneeSpeed = [simResultsT1.R.knee_speed.data; simResultsT2.R.knee_speed.data];
s1conc.RankleAngle = [simResultsT1.R.ankle_angle.data; simResultsT2.R.ankle_angle.data];
s1conc.RankleSpeed = [simResultsT1.R.ankle_speed.data; simResultsT2.R.ankle_speed.data];

s1conc.LhipAngle = [simResultsT1.L.hip_angle.data; simResultsT2.L.hip_angle.data];
s1conc.LhipSpeed = [simResultsT1.L.hip_speed.data; simResultsT2.L.hip_speed.data];
s1conc.LkneeAngle = [simResultsT1.L.knee_angle.data; simResultsT2.L.knee_angle.data];
s1conc.LkneeSpeed = [simResultsT1.L.knee_speed.data; simResultsT2.L.knee_speed.data];
s1conc.LankleAngle = [simResultsT1.L.ankle_angle.data; simResultsT2.L.ankle_angle.data];
s1conc.LankleSpeed = [simResultsT1.L.ankle_speed.data; simResultsT2.L.ankle_speed.data];

%Interpolate forces to match simulation results
gaitPeriod = length(subj1T1.forces.FP1(:,3))/subj1T1.framerate{1};
time = linspace(0,gaitPeriod,length(subj1T1.km.Rankle))';  

forceinterpT1 = interp1(time,subj1T1.forces.FP1(:,3),simResultsT1.R.ankle_angle.time);

gaitPeriod = length(subj1T2.forces.FP1(:,3))/subj1T2.framerate{1};
time = linspace(0,gaitPeriod,length(subj1T2.km.Rankle))';  

forceinterpT2 = interp1(time,subj1T2.forces.FP1(:,3),simResultsT2.R.ankle_angle.time);

s1conc.FP1 = [forceinterpT1; forceinterpT2];

%% Build regression variable (current method: Throw it all at the wall and see what sticks)

regressionVarSim = [s1conc.pelvis.px, s1conc.pelvis.vx, s1conc.pelvis.py, s1conc.pelvis.vy, s1conc.pelvis.pz, s1conc.pelvis.vz,...
    s1conc.pelvis.qx, s1conc.pelvis.qy, s1conc.pelvis.qz,s1conc.RhipAngle, s1conc.RhipSpeed, s1conc.RkneeAngle,...
    s1conc.RkneeSpeed, s1conc.RankleAngle, s1conc.RankleSpeed, s1conc.LhipAngle, s1conc.LhipSpeed,...
    s1conc.LkneeAngle, s1conc.LkneeSpeed, s1conc.LankleAngle, s1conc.LankleSpeed, s1conc.FP1];

% 
% RegressionVarExp = [subj1.km.Rhip(:,2),subj1.km.Rknee(:,2),subj1.km.PelvisAngle(:,1),...
%     subj1.km.PelvisAngle(:,2),subj1.km.PelvisAngle(:,3),subj1.km.PelvisCOMPos(:,1),...
%     subj1.km.PelvisCOMPos(:,3),subj1.forces.FP1(:,3)];
% 
% RegressionVarExp = [subj1.km.Rhip(:,2),subj1.km.Rknee(:,2),subj1.km.PelvisAngle(:,1),...
%     subj1.km.PelvisAngle(:,2),subj1.km.PelvisAngle(:,3),subj1.km.PelvisCOMPos(:,1),...
%     subj1.km.PelvisCOMPos(:,3)];

%% Assess model response for third trial
%Create regression variable for third trial
regressionVarT3 = [simResultsT3.pelvis.px.data, simResultsT3.pelvis.vx.data, simResultsT3.pelvis.py.data, simResultsT3.pelvis.vy.data,...
    simResultsT3.pelvis.pz.data, simResultsT3.pelvis.vz.data, simResultsT3.pelvis.qx.data, simResultsT3.pelvis.qy.data,...
    simResultsT3.pelvis.qz.data, simResultsT3.R.hip_angle.data, simResultsT3.R.hip_speed.data, simResultsT3.R.knee_angle.data,...
    simResultsT3.R.knee_speed.data, simResultsT3.R.ankle_angle.data, simResultsT3.R.ankle_speed.data,simResultsT3.L.hip_angle.data,...
    simResultsT3.L.hip_speed.data, simResultsT3.L.knee_angle.data, simResultsT3.L.knee_speed.data, simResultsT3.L.ankle_angle.data,...
    simResultsT3.L.ankle_speed.data];


% [trainedModel, validationRMSE] = trainRegressionModel1(RegressionVarSim);%Train model with new data

yfit = trainedModel6.predictFcn(regressionVarT3);

% figure;
% plot(subj1.forces.FP1(:,3)); hold on;
% plot(trainedModel.RegressionGP.Y);

gaitPeriod = length(subj1T3.forces.FP1(:,3))/subj1T3.framerate{1};
time = linspace(0,gaitPeriod,length(subj1T3.km.Rankle))';  

forceinterpT3 = interp1(time,subj1T3.forces.FP1(:,3),simResultsT3.R.ankle_angle.time);


figure;
plot(forceinterpT3); hold on;
plot(yfit);