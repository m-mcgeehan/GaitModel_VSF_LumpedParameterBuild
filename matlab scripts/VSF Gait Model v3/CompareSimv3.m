%% Compare forces
VSFGRFR_exp = downsample(sqrt(dynamic.forces.FP2(200:475,3).^2 + dynamic.forces.FP2(200:475,2).^2 + ...
    dynamic.forces.FP2(200:475,1).^2),2); 
VSFtime_exp = downsample(time(200:475),2);
VSFGRFR_exp = timeseries(VSFGRFR_exp,VSFtime_exp);

intactGRFR_exp = downsample(sqrt(dynamic.forces.FP1(160:350,3).^2 + dynamic.forces.FP1(160:350,2).^2 ...
    + dynamic.forces.FP1(160:350,1).^2),2);
intactTime_exp = downsample(time(160:350),2);
intactGRFR_exp = timeseries(intactGRFR_exp,intactTime_exp);

VSFGRFR_sim = SDOSimTest_Log.VSF_CF_resultant;
intactGRFR_sim = SDOSimTest_Log.intact_CF;

VSFGRFR_interp = interp1(VSFGRFR_sim.time, VSFGRFR_sim.data,VSFGRFR_exp.time, 'linear', 'extrap'); 
VSF_RMSE = sqrt(immse(VSFGRFR_exp.data,VSFGRFR_interp));
VSFGRFR_r2 = rsquared(VSFGRFR_interp,VSFGRFR_exp.data,2);

intactGRFR_interp = interp1(intactGRFR_sim.time, intactGRFR_sim.data, intactGRFR_exp.time,'linear','extrap');
intact_RMSE = sqrt(immse(intactGRFR_exp.data,intactGRFR_interp));
intactGRFR_r2 = rsquared(intactGRFR_interp,intactGRFR_exp.data,2);

%Impulse
VSFGRFR_exp_impulse = cumtrapz(VSFGRFR_exp.time,VSFGRFR_exp.data);
VSFGRFR_sim_impulse = cumtrapz(VSFGRFR_exp.time,VSFGRFR_interp);
VSF_impulse_RMSE = sqrt(immse(VSFGRFR_exp_impulse,VSFGRFR_sim_impulse));
VSFGRFR_impulse_r2 = rsquared(VSFGRFR_exp_impulse,VSFGRFR_sim_impulse,2);

intact_exp_impulse = cumtrapz(intactGRFR_exp.time,intactGRFR_exp.data);
intact_sim_impulse = cumtrapz(intactGRFR_exp.time,intactGRFR_interp);
intact_impulse_RMSE = sqrt(immse(intact_exp_impulse,intact_sim_impulse));
intact_impulse_r2 = rsquared(intact_exp_impulse,intact_sim_impulse,2);

figure;
subplot(2,2,1); plot(VSFGRFR_exp.time, VSFGRFR_interp,'k'); hold on; plot(VSFGRFR_exp,'r');...
    legend('simulation','experimental'); title(['VSF GRFR R^2 = ', num2str(VSFGRFR_r2), ' RMSE =', num2str(VSF_RMSE), 'N']);
subplot(2,2,2); plot(intactGRFR_exp.time, intactGRFR_interp,'k'); hold on; plot(intactGRFR_exp,'r'); legend('simulation','experimental');
    legend('simulation','experimental'); title(['Intact GRFR R^2 = ', num2str(intactGRFR_r2), ' RMSE =', num2str(intact_RMSE), 'N']);
subplot(2,2,3); plot(VSFGRFR_exp.time, VSFGRFR_sim_impulse,'k'); hold on; plot(VSFGRFR_exp.time,VSFGRFR_exp_impulse,'r');...
    legend('simulation','experimental'); title(['VSF Impulse R^2 = ', num2str(VSFGRFR_impulse_r2), ' RMSE =', num2str(VSF_impulse_RMSE), 'N/m/s']);
subplot(2,2,4); plot(intactGRFR_exp.time, intact_sim_impulse,'k'); hold on; plot(intactGRFR_exp.time,intact_exp_impulse,'r');...
    legend('simulation','experimental'); title(['Intact Impulse R^2 = ', num2str(intact_impulse_r2), ' RMSE =', num2str(intact_impulse_RMSE), 'N/m/s']);
sgtitle(strcat('Resultant GRF for: ', dynamic.trial));

% plot(VSFGRFR_exp.time, VSFGRFR_interp,'k'); hold on; plot(VSFGRFR_exp,'r'); ...
%     hold on; plot(intactGRFR_exp.time, intactGRFR_interp,'k'); hold on; plot(intactGRFR_exp,'r');



%% Kinematics (sagittal plane)
idx = 450;

RAnkleSim = interp1(SDOSimTest_Log.tout, SDOSimTest_Log.measR.ankle_angle.data,time, 'linear', 'extrap');
RAnkle_RMSE = sqrt(immse(ankle_motionRY(1:idx,2),RAnkleSim(1:idx)));
RAnkle_r2 = rsquared(RAnkleSim(1:idx), ankle_motionRY(1:idx,2),2);
% plot(RAnkleSim(1:idx),'k'); hold on; plot(rad2deg(ankle_motionRY(1:idx,2)),'r');

LAnkleSim = interp1(SDOSimTest_Log.tout, SDOSimTest_Log.measL.ankle_angle.data,time, 'linear', 'extrap');
LAnkle_RMSE = sqrt(immse(rad2deg(ankleLAngY(1:idx,2)),LAnkleSim(1:idx)));
LAnkle_r2 = rsquared(LAnkleSim(1:idx), rad2deg(ankleLAngY(1:idx,2)),2);
% plot(LAnkleSim(1:idx),'k'); hold on; plot(rad2deg(ankleLAngY(1:idx,2)),'r');

RKneeSim = interp1(SDOSimTest_Log.tout, SDOSimTest_Log.measR.knee_angle.data,time, 'linear', 'extrap');
RKnee_RMSE = sqrt(immse(rad2deg(kneeRAngY(1:idx,2)),RKneeSim(1:idx)));
RKnee_r2 = rsquared(RKneeSim(1:idx), rad2deg(kneeRAngY(1:idx,2)),2);
% plot(RKneeSim(1:idx),'k'); hold on; plot(rad2deg(kneeRAngY(1:idx,2)),'r');

LKneeSim = interp1(SDOSimTest_Log.tout, SDOSimTest_Log.measL.knee_angle.data,time, 'linear', 'extrap');
LKnee_RMSE = sqrt(immse(rad2deg(kneeLAngY(1:idx,2)),LKneeSim(1:idx)));
LKnee_r2 = rsquared(LKneeSim(1:idx), rad2deg(kneeLAngY(1:idx,2)),2);
% plot(LKneeSim(1:idx),'k'); hold on; plot(rad2deg(kneeLAngY(1:idx,2)),'r');

RHipSim = interp1(SDOSimTest_Log.tout, SDOSimTest_Log.measR.hip_angle.data,time, 'linear', 'extrap');
RHip_RMSE = sqrt(immse(rad2deg(hipRAngY(1:idx,2)),RHipSim(1:idx)));
RHip_r2 = rsquared(RHipSim(1:idx), rad2deg(hipRAngY(1:idx,2)),2);
% plot(RHipSim(1:idx),'k'); hold on; plot(rad2deg(hipRAngY(1:idx,2)),'r');

LHipSim = interp1(SDOSimTest_Log.tout, SDOSimTest_Log.measL.hip_angle.data,time, 'linear', 'extrap');
LHip_RMSE = sqrt(immse(rad2deg(hipLAngY(1:idx,2)),LHipSim(1:idx)));
LHip_r2 = rsquared(LHipSim(1:idx), rad2deg(hipLAngY(1:idx,2)),2);
% plot(LHipSim(1:idx),'k'); hold on; plot(rad2deg(hipLAngY(1:idx,2)),'r');

pelvisXSim = interp1(SDOSimTest_Log.tout, 100*SDOSimTest_Log.measPelvis.px.data,time, 'linear', 'extrap');
pelvisX_RMSE = sqrt(immse(pelvisPosX(1:idx,2),pelvisXSim(1:idx)));
pelvisX_r2 = rsquared(pelvisPosX(1:idx,2),pelvisXSim(1:idx));

% 
% pelvisYSim = interp1(SDOSimTest_Log.tout, 100*SDOSimTest_Log.measPelvis.py.data,time, 'linear', 'extrap');
% pelvisY_RMSE = sqrt(immse(pelvisPosY(1:idx,2),pelvisYSim(1:idx)));
% pelvisY_r2 = rsquared(pelvisPosY(1:idx,2),pelvisYSim(1:idx));


figure;
subplot(4,2,1); plot(LHipSim(1:idx),'k'); hold on; plot(rad2deg(hipLAngY(1:idx,2)),'r');
    title('LHip');legend('simulation','experimental'); title(['R^2 = ', num2str(LHip_r2), ' RMSE = ', num2str(LHip_RMSE), 'deg']); ylabel('HipL Angle (deg)');
subplot(4,2,2); plot(RHipSim(1:idx),'k'); hold on; plot(rad2deg(hipRAngY(1:idx,2)),'r');
    title('RHip');legend('simulation','experimental'); title(['R^2 = ', num2str(RHip_r2), ' RMSE = ', num2str(RHip_RMSE), 'deg']); ylabel('HipR Angle (deg)');
subplot(4,2,3); plot(LKneeSim(1:idx),'k'); hold on; plot(rad2deg(kneeLAngY(1:idx,2)),'r');
    title('LKnee');legend('simulation','experimental'); title(['R^2 = ', num2str(LKnee_r2), ' RMSE = ', num2str(LKnee_RMSE), 'deg']); ylabel('KneeL Angle (deg)');
subplot(4,2,4); plot(RKneeSim(1:idx),'k'); hold on; plot(rad2deg(kneeRAngY(1:idx,2)),'r');
    title('RKnee');legend('simulation','experimental'); title(['R^2 = ', num2str(RKnee_r2), ' RMSE = ', num2str(RKnee_RMSE), 'deg']); ylabel('KneeR Angle (deg)');
subplot(4,2,5); plot(LAnkleSim(1:idx),'k'); hold on; plot(rad2deg(ankleLAngY(1:idx,2)),'r');
    title('RAnkle');legend('simulation','experimental'); title(['R^2 = ', num2str(LAnkle_r2), ' RMSE = ', num2str(LAnkle_RMSE), 'deg']);ylabel('AnkleL Angle (deg)');
subplot(4,2,6); plot(RAnkleSim(1:idx),'k'); hold on; plot(rad2deg(ankle_motionRY(1:idx,2)),'r');
    title('RAnkle');legend('simulation','experimental'); title(['R^2 = ', num2str(RAnkle_r2), ' RMSE = ', num2str(RAnkle_RMSE), 'deg']); ylabel('AnkleR Angle (deg)');
subplot(4,2,7); plot(pelvisXSim(1:idx),'k'); hold on; plot(pelvisPosX(1:idx,2),'r');
    title('PelvisPosX');legend('simulation','experimental'); title(['R^2 = ', num2str(pelvisX_r2), ' RMSE = ', num2str(pelvisX_RMSE), 'cm']);ylabel('Pelvis PosX (cm)');
sgtitle(strcat('Sagittal plane kinematics for: ', dynamic.trial));



% 
% fc = 6;%6 hz low pass cut off
% Wn_fp = 2*fc/fs; 
% [b,a] = butter(order,Wn_fp);
% 
% PelvisPosFilt = filtfilt(b,a,measPelvis.px.data);
% measPelvisPos.Pos.data = PelvisPosFilt; clear PelvisPosFilt
% 
% RHipAngleFilt = filtfilt(b,a,measR.hip_angle.data);
% measR.hip_angle.data = RHipAngleFilt; clear RHipAngleFilt
% 
% LHipAngleFilt = filtfilt(b,a,measL.hip_angle.data);
% measL.hip_angle.data = LHipAngleFilt; clear LHipAngleFilt
% 
% RKneeAngleFilt = filtfilt(b,a,measR.knee_angle.data);
% measR.knee_angle.data = RKneeAngleFilt; clear RKneeAngleFilt
% 
% LKneeAngleFilt = filtfilt(b,a,measL.knee_angle.data);
% measL.knee_angle.data = LKneeAngleFilt; clear LKneeAngleFilt
% 
% RAnkleAngleFilt = filtfilt(b,a,measR.ankle_angle.data);
% measR.ankle_angle.data = RAnkleAngleFilt; clear RAnkleAngleFilt
% 
% LAnkleAngleFilt = filtfilt(b,a,measL.ankle_angle.data);
% measL.ankle_angle.data = LAnkleAngleFilt; clear LAnkleAngleFilt
% 
% % PelvisAngleFilt = filtfilt(b,a,PelvisAngle);
% % PelvisCOMPosFilt = filtfilt(b,a,PelvisCOMPos);
% % 
% 
% time1 = linspace(0,gaitPeriod,length(dynamic.km.Rhip(:,3)));
% figure;  
% subplot(3,3,1);plot(timeSim,(measR.hip_angle.data),'r');hold on;plot(time1,rad2deg(hipRAngY(:,2)),'k');
% title('Right Hip angle')
% subplot(3,3,2);plot(timeSim,(measL.hip_angle.data),'r');hold on;plot(time1,rad2deg(hipLAngY(:,2)),'k');
% title('Left Hip angle')
% subplot(3,3,3);plot(timeSim,(measR.knee_angle.data),'r');hold on;plot(time1,rad2deg(knee_motionRY(:,2)),'k');
% title('Right Knee angle')
% subplot(3,3,4);plot(timeSim,(measL.knee_angle.data),'r');hold on;plot(time1,rad2deg(knee_motionLY(:,2)),'k');
% title('Left Knee angle')
% subplot(3,3,5);plot(timeSim,(measR.ankle_angle.data),'r');hold on;plot(time1,rad2deg(ankle_motionRY(:,2)),'k');
% title('Right Ankle angle')
% subplot(3,3,6);plot(timeSim,(measL.ankle_angle.data),'r');hold on;plot(time1,rad2deg(ankle_motionLY(:,2)),'k');
% title('Left Ankle angle')
% subplot(3,3,7);plot(timeSim,(measPelvis.px.data),'r');hold on;plot(time1,pelvisPosX(:,2),'k');
% title('Pelvis Position')
% sgtitle('Simulation vs Experimental kinematics (red = sim, black = exp)');
% 
% 
% figure;
% subplot(2,3,1);plot(time1,rad2deg(pelvisAngX(:,2)),'k');hold on;plot(timeSim,rad2deg(measPelvis.qx.data),'r');
% title('Pelvis angleX')
% subplot(2,3,2);plot(time1,rad2deg(pelvisAngY(:,2)),'k');hold on;plot(timeSim,rad2deg(measPelvis.qy.data),'r');
% title('Pelvis angleY')
% subplot(2,3,3);plot(time1,rad2deg(pelvisAngZ(:,2)),'k');hold on;plot(timeSim,rad2deg(measPelvis.qz.data),'r');
% title('Pelvis angleZ')
% subplot(2,3,4);plot(time1,pelvisPosX(:,2),'k');hold on;plot(timeSim,measPelvis.px.data,'r');
% title('Pelvis PosX')
% subplot(2,3,5);plot(time1,pelvisPosY(:,2),'k');hold on;plot(timeSim,measPelvis.py.data,'r');
% title('Pelvis PosY')
% subplot(2,3,6);plot(time1,pelvisPosZ(:,2),'k');hold on;plot(timeSim,measPelvis.pz.data,'r');
% title('Pelvis PosZ')
% sgtitle('Simulation vs Experimental Pelvis kinematics (red = sim, black = exp)');
% 
% 
% %% 
% a = linspace(0,gaitPeriod,length(dynamic.forces.FP1));
% 
% FP1TS = timeseries(dynamic.forces.FP1,a);
% FP2TS = timeseries(dynamic.forces.FP2,a);
% 
% order = 2;
% fs = length(VSF_CF.time)/gaitPeriod; 
% fc = 15;%50 hz cut off
% Wn_fp = 2*fc/fs; 
% [bfp,afp] = butter(order,Wn_fp); 
% simForceVSF = filtfilt(bfp,afp,VSF_CF.data); %symmetric lowpass filter of normal force
% simForceIntact = filtfilt(bfp,afp,intact_CF.data); %symmetric lowpass filter of normal force
% timeSim = linspace(0,gaitPeriod,length(measR.normal_force.time));
% 
% 
% plot(FP1TS.time,FP1TS.data(:,3)); hold on; plot(FP2TS.time,FP2TS.data(:,3)); hold on;...
%     plot(VSF_CF.time,simForceVSF); hold on; plot(intact_CF.time,simForceIntact); hold on; xlim([.85 2.25]);

% %Force
% order = 2;
% fs = length(VSF_CF.time)/gaitPeriod; 
% fc = 20;%50 hz cut off
% Wn_fp = 2*fc/fs; 
% [bfp,afp] = butter(order,Wn_fp); 
% simForceR = filtfilt(bfp,afp,VSF_CF.data); %symmetric lowpass filter of normal force
% simForceL = filtfilt(bfp,afp,intact_CF.data); %symmetric lowpass filter of normal force
% timeSim = linspace(0,gaitPeriod,length(VSF_CF.time));
% 
% time1 = linspace(0,gaitPeriod,length(dynamic.forces.FP1(:,3)));
% figure; hold on
% plot(timeSim,simForceR); hold on
% plot(timeSim,simForceL); hold on
% plot(time1,dynamic.forces.FP1(:,3)); hold on
% plot(time1,dynamic.forces.FP2(:,3)); hold on
% xlim([.85 2.25]);

% 
% figure;  
% subplot(3,2,1);plot(timeSim,rad2deg(measR.hip_angle.data),'r');hold on;plot(time,dynamic.km.Rhip(:,2),'k');
% title('Right Hip angle')
% subplot(3,2,2);plot(timeSim,rad2deg(measL.hip_angle.data),'r');hold on;plot(time,dynamic.km.Lhip(:,2),'k');
% title('Left Hip angle')
% subplot(3,2,3);plot(timeSim,rad2deg(measR.knee_angle.data),'r');hold on;plot(time,dynamic.km.Rknee(:,2),'k');
% title('Right Knee angle')
% subplot(3,2,4);plot(timeSim,rad2deg(measL.knee_angle.data),'r');hold on;plot(time,dynamic.km.Lknee(:,2),'k');
% title('Left Knee angle')
% subplot(3,2,5);plot(timeSim,rad2deg(measR.ankle_angle.data),'r');hold on;plot(time,dynamic.km.Rankle(:,2),'k');
% title('Right Ankle angle')
% subplot(3,2,6);plot(timeSim,rad2deg(measL.ankle_angle.data),'r');hold on;plot(time,dynamic.km.Lankle(:,2),'k');
% title('Left Ankle angle')
% sgtitle('Simulation vs Experimental kinematics (red = sim, black = exp)');