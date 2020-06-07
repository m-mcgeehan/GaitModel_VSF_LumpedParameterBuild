%% Compare forces
saveData = 1; %SaveData? 1 = yes, 0 = no
fs = 200; %hz
buffer = fs*.25; %number of samples to pad on either side of stance onset/offset
time = linspace(0,(length(dynamic.forces.FP2)/fs),length(dynamic.forces.FP2));
phaseShift = .02; %phase shift (s) from analog low-pass filter at 40 Hz in simulation

%%
idx1 = find(dynamic.forces.FP2(:,3) > 1,1,'first');
idx2 = find(dynamic.forces.FP2(:,3)> 1,1,'last');

VSFGRFR_exp = timeseries((sqrt(dynamic.forces.FP2(:,3).^2 + dynamic.forces.FP2(:,2).^2 + ...
    dynamic.forces.FP2(:,1).^2)),time+phaseShift);

% enddif = (length(VSFGRFR_exp.data(idx2:end)))+buffer;
% 
% if  (length(VSFGRFR_exp.data(idx1-buffer:idx2+buffer))) < (length(VSFGRFR_exp.data))
%     
% else
%     VSFGRFR_exp = extendts(VSFGRFR_exp,(length(VSFGRFR_exp.data)+50),0);
% end

% if (length(VSFGRFR_exp.data(idx1-buffer:idx2+buffer))) > (length(VSFGRFR_exp.data)) %Check length requirements
%     warning('Not enough experimental data for .25 s pad')
% else
% end


a = length(VSFGRFR_exp.data(idx2:end));

if a < buffer %Check length requirements
    warning('Not enough experimental data for .25 s pad...padding time series data with zeros')
    pad = timeseries((linspace(0,0,buffer)),(linspace(VSFGRFR_exp.time(end),VSFGRFR_exp.time(end)+.25,buffer)));
    VSFGRFR_exp = append(VSFGRFR_exp,pad);
else
end

% VSFGRFR_expShift2 = timeseries(downsample(VSFGRFR_exp.data(idx1-buffer:idx2+buffer),2),...
%      downsample((VSFGRFR_exp.time(idx1-buffer:idx2+buffer)+phaseShift),2));
   
VSFGRFR_expShift = timeseries(downsample(VSFGRFR_exp.data(idx1-buffer:idx2+buffer),2),...
    downsample(VSFGRFR_exp.time(idx1-buffer:idx2+buffer),2));


idx1 = find(dynamic.forces.FP1(:,3) > 1,1,'first');
idx2 = find(dynamic.forces.FP1(:,3)> 1,1,'last');
intactGRFR_exp = timeseries((sqrt(dynamic.forces.FP1(:,3).^2 + dynamic.forces.FP1(:,2).^2 + ...
    dynamic.forces.FP1(:,1).^2)),time);

% if  (length(intactGRFR_exp.data(idx1-buffer:idx2+buffer))) < (length(intactGRFR_exp.data))
%     
% else
%     VSFGRFR_exp = extendts(intactGRFR_exp,(length(intactGRFR_exp.data)+50),0);
% end

intactGRFR_expShift = timeseries(downsample(intactGRFR_exp.data(idx1-buffer:idx2+buffer),2),...
    downsample(intactGRFR_exp.time(idx1-buffer:idx2+buffer),2));


VSFGRFR_sim = SDOSimTest_Log.VSF_CF_resultant;
intactGRFR_sim = SDOSimTest_Log.intact_CF;

VSFGRFR_interp = interp1(VSFGRFR_sim.time, VSFGRFR_sim.data,VSFGRFR_expShift.time, 'spline', 'extrap'); 
VSF_RMSE = sqrt(immse(VSFGRFR_expShift.data,VSFGRFR_interp));
VSFGRFR_r2 = rsquared(VSFGRFR_interp,VSFGRFR_expShift.data,2);

intactGRFR_interp = interp1(intactGRFR_sim.time, intactGRFR_sim.data, intactGRFR_expShift.time,'spline','extrap');
intact_RMSE = sqrt(immse(intactGRFR_expShift.data,intactGRFR_interp));
intactGRFR_r2 = rsquared(intactGRFR_interp,intactGRFR_expShift.data,2);



%Impulse
VSFGRFR_exp_impulse = cumtrapz(VSFGRFR_expShift.time,VSFGRFR_expShift.data);
VSFGRFR_sim_impulse = cumtrapz(VSFGRFR_expShift.time,VSFGRFR_interp);
VSF_impulse_RMSE = sqrt(immse(VSFGRFR_exp_impulse,VSFGRFR_sim_impulse));
VSFGRFR_impulse_r2 = rsquared(VSFGRFR_exp_impulse,VSFGRFR_sim_impulse,2);

intact_exp_impulse = cumtrapz(intactGRFR_expShift.time,intactGRFR_expShift.data);
intact_sim_impulse = cumtrapz(intactGRFR_expShift.time,intactGRFR_interp);
intact_impulse_RMSE = sqrt(immse(intact_exp_impulse,intact_sim_impulse));
intact_impulse_r2 = rsquared(intact_exp_impulse,intact_sim_impulse,2);

figure;
subplot(2,2,1); plot(VSFGRFR_expShift.time, VSFGRFR_interp,'k'); hold on; plot(VSFGRFR_expShift,'r');...
    legend('simulation','experimental'); title(['VSF GRFR R^2 = ', num2str(VSFGRFR_r2), ' RMSE =', num2str(VSF_RMSE), 'N']);
subplot(2,2,2); plot(intactGRFR_expShift.time, intactGRFR_interp,'k'); hold on; plot(intactGRFR_expShift,'r'); legend('simulation','experimental');
    legend('simulation','experimental'); title(['Intact GRFR R^2 = ', num2str(intactGRFR_r2), ' RMSE =', num2str(intact_RMSE), 'N']);
subplot(2,2,3); plot(VSFGRFR_expShift.time, VSFGRFR_sim_impulse,'k'); hold on; plot(VSFGRFR_expShift.time,VSFGRFR_exp_impulse,'r');...
    legend('simulation','experimental'); title(['VSF Impulse R^2 = ', num2str(VSFGRFR_impulse_r2), ' RMSE =', num2str(VSF_impulse_RMSE), 'N/dots']);
subplot(2,2,4); plot(intactGRFR_expShift.time, intact_sim_impulse,'k'); hold on; plot(intactGRFR_expShift.time,intact_exp_impulse,'r');...
    legend('simulation','experimental'); title(['Intact Impulse R^2 = ', num2str(intact_impulse_r2), ' RMSE =', num2str(intact_impulse_RMSE), 'N/dots']);
sgtitle(strcat('Resultant GRF for: ', dynamic.trial));


    

%% Kinematics (sagittal plane)
idx = 450;
time = linspace(0,(length(dynamic.forces.FP2)/fs),length(dynamic.forces.FP2));
RAnkleSim = interp1(SDOSimTest_Log.tout, SDOSimTest_Log.measR.ankle_angle.data,time, 'spline', 'extrap');
RAnkle_RMSE = sqrt(immse(ankle_motionRY(1:idx,2),RAnkleSim(1:idx)'));
RAnkle_r2 = rsquared(RAnkleSim(1:idx)', ankle_motionRY(1:idx,2),2);
% plot(RAnkleSim(1:idx),'k'); hold on; plot(rad2deg(ankle_motionRY(1:idx,2)),'r');

LAnkleSim = interp1(SDOSimTest_Log.tout, SDOSimTest_Log.measL.ankle_angle.data,time, 'spline', 'extrap');
LAnkle_RMSE = sqrt(immse(rad2deg(ankleLAngY(1:idx,2)),LAnkleSim(1:idx)'));
LAnkle_r2 = rsquared(LAnkleSim(1:idx)', rad2deg(ankleLAngY(1:idx,2)),2);
% plot(LAnkleSim(1:idx),'k'); hold on; plot(rad2deg(ankleLAngY(1:idx,2)),'r');

RKneeSim = interp1(SDOSimTest_Log.tout, SDOSimTest_Log.measR.knee_angle.data,time, 'spline', 'extrap');
RKnee_RMSE = sqrt(immse(rad2deg(kneeRAngY(1:idx,2)),RKneeSim(1:idx)'));
RKnee_r2 = rsquared(RKneeSim(1:idx)', rad2deg(kneeRAngY(1:idx,2)),2);
% plot(RKneeSim(1:idx),'k'); hold on; plot(rad2deg(kneeRAngY(1:idx,2)),'r');

LKneeSim = interp1(SDOSimTest_Log.tout, SDOSimTest_Log.measL.knee_angle.data,time, 'spline', 'extrap');
LKnee_RMSE = sqrt(immse(rad2deg(kneeLAngY(1:idx,2)),LKneeSim(1:idx)'));
LKnee_r2 = rsquared(LKneeSim(1:idx)', rad2deg(kneeLAngY(1:idx,2)),2);
% plot(LKneeSim(1:idx),'k'); hold on; plot(rad2deg(kneeLAngY(1:idx,2)),'r');

RHipSim = interp1(SDOSimTest_Log.tout, SDOSimTest_Log.measR.hip_angle.data,time, 'spline', 'extrap');
RHip_RMSE = sqrt(immse(rad2deg(hipRAngY(1:idx,2)),RHipSim(1:idx)'));
RHip_r2 = rsquared(RHipSim(1:idx)', rad2deg(hipRAngY(1:idx,2)),2);
% plot(RHipSim(1:idx),'k'); hold on; plot(rad2deg(hipRAngY(1:idx,2)),'r');

LHipSim = interp1(SDOSimTest_Log.tout, SDOSimTest_Log.measL.hip_angle.data,time, 'spline', 'extrap');
LHip_RMSE = sqrt(immse(rad2deg(hipLAngY(1:idx,2)),LHipSim(1:idx)'));
LHip_r2 = rsquared(LHipSim(1:idx)', rad2deg(hipLAngY(1:idx,2)),2);
% plot(LHipSim(1:idx),'k'); hold on; plot(rad2deg(hipLAngY(1:idx,2)),'r');

pelvisXSim = interp1(SDOSimTest_Log.tout, 100*SDOSimTest_Log.measPelvis.px.data,time, 'spline', 'extrap');
pelvisX_RMSE = sqrt(immse(pelvisPosX(1:idx,2),pelvisXSim(1:idx)'));
pelvisX_r2 = rsquared(pelvisPosX(1:idx,2),pelvisXSim(1:idx)');

% 
% pelvisYSim = interp1(SDOSimTest_Log.tout, 100*SDOSimTest_Log.measPelvis.py.data,time, 'spline', 'extrap');
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
%% Calc CoP
% calcCOP

%% Save data?
if isequal(saveData,1)
    
    struct.dynamic = dynamic;
    struct.SDOSimTest_Log = SDOSimTest_Log;
    struct.compare.VSFGRFR_expShift = VSFGRFR_expShift;
    struct.compare.intactGRFR_expShift =intactGRFR_expShift;
    struct.compare.VSFGRFR_interp = VSFGRFR_interp;
    struct.compare.VSF_RMSE = VSF_RMSE;
    struct.compare.VSFGRFR_r2 = VSFGRFR_r2;
    struct.compare.intact_RMSE = intact_RMSE;
    struct.compare.intactGRFR_r2 = intactGRFR_r2;
    struct.compare.VSFGRFR_exp_impulse = VSFGRFR_exp_impulse;
    struct.compare.VSFGRFR_sim_impulse = VSFGRFR_sim_impulse;
    struct.compare.VSF_impulse_RMSE = VSF_impulse_RMSE;
    struct.compare.VSFGRFR_impulse_r2 = VSFGRFR_impulse_r2;
    struct.compare.intact_exp_impulse = intact_exp_impulse;
    struct.compare.intact_sim_impulse = intact_sim_impulse;
    struct.compare.intact_impulse_RMSE = intact_impulse_RMSE;
    struct.compare.intact_impulse_r2 = intact_impulse_r2;
    struct.compare.RAnkle_RMSE = RAnkle_RMSE;
    struct.compare.RAnkle_r2 = RAnkle_r2;
    struct.compare.LAnkle_RMSE = LAnkle_RMSE;
    struct.compare.LAnkle_r2 = LAnkle_r2;
    struct.compare.RKnee_RMSE = RKnee_RMSE;
    struct.compare.RKnee_r2 = RKnee_r2;
    struct.compare.LKnee_RMSE = LKnee_RMSE;
    struct.compare.LKnee_r2 = LKnee_r2;
    struct.compare.RHip_RMSE = RHip_RMSE;
    struct.compare.RHip_r2 = RHip_r2;
    struct.compare.LHip_RMSE = LHip_RMSE;
    struct.compare.LHip_r2 = LHip_r2;
    struct.compare.pelvisX_RMSE = pelvisX_RMSE;
    struct.compare.pelvisX_r2 = pelvisX_r2;
    struct.compare.COP.copX101 = copX101;
    struct.compare.COP.VSFCOP_exp101 = VSFCOP_exp101;
    struct.compare.COP.copX = copX;
    struct.compare.COP.VSFCOP_exp = VSFCOP_exp;
    struct.compare.COP.COPRMSE = copX_RMSE;
    struct.compare.COP.COP_r2 = copX_r2;
    uisave('struct');
    uisave();
else
end

%%
% idx1 = find(struct.compare.VSFGRFR_expShift.data > 1,1,'first');
% idx2 = find(struct.compare.VSFGRFR_expShift.data > 1,1,'last');
% stanceTime_VSFexp = struct.compare.VSFGRFR_expShift.time(idx2) - struct.compare.VSFGRFR_expShift.time(idx1);
% 
% plot(struct.compare.VSFGRFR_expShift); hold on; plot(struct.compare.VSFGRFR_expShift.time,struct.compare.VSFGRFR_interp); 
% 
% idx1 = find(struct.compare.VSFGRFR_interp > 1,1,'first');
% 
% pause('on');
% 
% %%
% pause('off');
% a = 2.274
% stanceTime_VSFsim = a - struct.compare.VSFGRFR_expShift.time(idx1);
% 
% 
