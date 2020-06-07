%% low stiffness
%Interpolate to 101 data points (stance -.25 --> stance +.25 s)
fs = 200; %hz
buffer = fs*.25; %number of samples to pad on either side of stance onset/offset


numTrials = 3; %Number of trials to generate curves for
subjNumber = 'VSF004';
% 
% for i = 1:numTrials
%     eval(['num =',num2str(i),';']);
%     
%     
%     
% end
% 
% 
% for i = 1:length(markerNames)
%     marker = markerNames{i};
%     eval(['markerData =',marker,';']);
%     markerDataIdx = markerData{1}(:,1:3)';
%     markerDataFilt = filtfilt(b,a,markerDataIdx');
%     markerDataIdx = markerDataFilt';
%     for k = 1:length(markerDataIdx(1,:))
%         markerDataRot(k,:) = (x*(markerDataIdx(:,k))*100);
%     end
%     markerSet.(marker) = markerDataRot;
%     markerSetTS.(marker) = timeseries(markerDataRot,time);
% end
% 


% GRFR
VSFs1t1_exp = VSF004_CAS1T1.compare.VSFGRFR_expShift;
VSFs1t2_exp = VSF004_CAS1T2.compare.VSFGRFR_expShift;
VSFs1t3_exp = VSF004_CAS1T3.compare.VSFGRFR_expShift;
VSFs1t1_sim = timeseries(VSF004_CAS1T1.compare.VSFGRFR_interp,VSF004_CAS1T1.compare.VSFGRFR_expShift.time);
VSFs1t2_sim = timeseries(VSF004_CAS1T2.compare.VSFGRFR_interp,VSF004_CAS1T2.compare.VSFGRFR_expShift.time);
VSFs1t3_sim = timeseries(VSF004_CAS1T3.compare.VSFGRFR_interp,VSF004_CAS1T3.compare.VSFGRFR_expShift.time);

VSFs2t1_exp = VSF004_CCS2T1.compare.VSFGRFR_expShift;
VSFs2t2_exp = VSF004_CCS2T2.compare.VSFGRFR_expShift;
VSFs2t3_exp = VSF004_CCS2T3.compare.VSFGRFR_expShift;
VSFs2t1_sim = timeseries(VSF004_CCS2T1.compare.VSFGRFR_interp,VSF004_CCS2T1.compare.VSFGRFR_expShift.time);
VSFs2t2_sim = timeseries(VSF004_CCS2T2.compare.VSFGRFR_interp,VSF004_CCS2T2.compare.VSFGRFR_expShift.time);
VSFs2t3_sim = timeseries(VSF004_CCS2T3.compare.VSFGRFR_interp,VSF004_CCS2T3.compare.VSFGRFR_expShift.time);

VSFs3t1_exp = VSF004_CBS3T1.compare.VSFGRFR_expShift;
VSFs3t2_exp = VSF004_CBS3T2.compare.VSFGRFR_expShift;
VSFs3t3_exp = VSF004_CBS3T3.compare.VSFGRFR_expShift;
VSFs3t1_sim = timeseries(VSF004_CBS3T1.compare.VSFGRFR_interp,VSF004_CBS3T1.compare.VSFGRFR_expShift.time);
VSFs3t2_sim = timeseries(VSF004_CBS3T2.compare.VSFGRFR_interp,VSF004_CBS3T2.compare.VSFGRFR_expShift.time);
VSFs3t3_sim = timeseries(VSF004_CBS3T3.compare.VSFGRFR_interp,VSF004_CBS3T3.compare.VSFGRFR_expShift.time);


%Impulse
VSFs1t1_impulse_exp = timeseries(VSF004_CAS1T1.compare.VSFGRFR_exp_impulse,VSF004_CAS1T1.compare.VSFGRFR_expShift.time);
VSFs1t2_impulse_exp = timeseries(VSF004_CAS1T2.compare.VSFGRFR_exp_impulse,VSF004_CAS1T2.compare.VSFGRFR_expShift.time);
VSFs1t3_impulse_exp = timeseries(VSF004_CAS1T3.compare.VSFGRFR_exp_impulse,VSF004_CAS1T3.compare.VSFGRFR_expShift.time);
VSFs1t1_impulse_sim = timeseries(VSF004_CAS1T1.compare.VSFGRFR_sim_impulse,VSF004_CAS1T1.compare.VSFGRFR_expShift.time);
VSFs1t2_impulse_sim = timeseries(VSF004_CAS1T2.compare.VSFGRFR_sim_impulse,VSF004_CAS1T2.compare.VSFGRFR_expShift.time);
VSFs1t3_impulse_sim = timeseries(VSF004_CAS1T3.compare.VSFGRFR_sim_impulse,VSF004_CAS1T3.compare.VSFGRFR_expShift.time);

VSFs2t1_impulse_exp = timeseries(VSF004_CCS2T1.compare.VSFGRFR_exp_impulse,VSF004_CCS2T1.compare.VSFGRFR_expShift.time);
VSFs2t2_impulse_exp = timeseries(VSF004_CCS2T2.compare.VSFGRFR_exp_impulse,VSF004_CCS2T2.compare.VSFGRFR_expShift.time);
VSFs2t3_impulse_exp = timeseries(VSF004_CCS2T3.compare.VSFGRFR_exp_impulse,VSF004_CCS2T3.compare.VSFGRFR_expShift.time);
VSFs2t1_impulse_sim = timeseries(VSF004_CCS2T1.compare.VSFGRFR_sim_impulse,VSF004_CCS2T1.compare.VSFGRFR_expShift.time);
VSFs2t2_impulse_sim = timeseries(VSF004_CCS2T2.compare.VSFGRFR_sim_impulse,VSF004_CCS2T2.compare.VSFGRFR_expShift.time);
VSFs2t3_impulse_sim = timeseries(VSF004_CCS2T3.compare.VSFGRFR_sim_impulse,VSF004_CCS2T3.compare.VSFGRFR_expShift.time);

VSFs3t1_impulse_exp = timeseries(VSF004_CBS3T1.compare.VSFGRFR_exp_impulse,VSF004_CBS3T1.compare.VSFGRFR_expShift.time);
VSFs3t2_impulse_exp = timeseries(VSF004_CBS3T2.compare.VSFGRFR_exp_impulse,VSF004_CBS3T2.compare.VSFGRFR_expShift.time);
VSFs3t3_impulse_exp = timeseries(VSF004_CBS3T3.compare.VSFGRFR_exp_impulse,VSF004_CBS3T3.compare.VSFGRFR_expShift.time);
VSFs3t1_impulse_sim = timeseries(VSF004_CBS3T1.compare.VSFGRFR_sim_impulse,VSF004_CBS3T1.compare.VSFGRFR_expShift.time);
VSFs3t2_impulse_sim = timeseries(VSF004_CBS3T2.compare.VSFGRFR_sim_impulse,VSF004_CBS3T2.compare.VSFGRFR_expShift.time);
VSFs3t3_impulse_sim = timeseries(VSF004_CBS3T3.compare.VSFGRFR_sim_impulse,VSF004_CBS3T3.compare.VSFGRFR_expShift.time);


trialNames = {'VSFt1_exp', 'VSFt2_exp', 'VSFt3_exp', 'VSFt1_sim', 'VSFt2_sim', 'VSFt3_sim'};

%% COP
VSFs1t1_COP_sim = VSF004_CAS1T1.compare.COP.copX101;
VSFs1t1_COP_exp = VSF004_CAS1T1.compare.COP.VSFCOP_exp101;
VSFs1t2_COP_sim = VSF004_CAS1T2.compare.COP.copX101;
VSFs1t2_COP_exp = VSF004_CAS1T2.compare.COP.VSFCOP_exp101;
VSFs1t3_COP_sim = VSF004_CAS1T3.compare.COP.copX101;
VSFs1t3_COP_exp = VSF004_CAS1T3.compare.COP.VSFCOP_exp101;

VSFs2t1_COP_sim = VSF004_CCS2T1.compare.COP.copX101;
VSFs2t1_COP_exp = VSF004_CCS2T1.compare.COP.VSFCOP_exp101;
VSFs2t2_COP_sim = VSF004_CCS2T2.compare.COP.copX101;
VSFs2t2_COP_exp = VSF004_CCS2T2.compare.COP.VSFCOP_exp101;
VSFs2t3_COP_sim = VSF004_CCS2T3.compare.COP.copX101;
VSFs2t3_COP_exp = VSF004_CCS2T3.compare.COP.VSFCOP_exp101;

VSFs3t1_COP_sim = VSF004_CBS3T1.compare.COP.copX101;
VSFs3t1_COP_exp = VSF004_CBS3T1.compare.COP.VSFCOP_exp101;
VSFs3t2_COP_sim = VSF004_CBS3T2.compare.COP.copX101;
VSFs3t2_COP_exp = VSF004_CBS3T2.compare.COP.VSFCOP_exp101;
VSFs3t3_COP_sim = VSF004_CBS3T3.compare.COP.copX101;
VSFs3t3_COP_exp = VSF004_CBS3T3.compare.COP.VSFCOP_exp101;
%% Refit Experimental data to 101 data points

VSFs1t1_exp = timeseries(imresize(VSFs1t1_exp.data,[101 1]),imresize(VSFs1t1_exp.time,[101 1])); %resample time series to 101 data points
VSFs1t1_impulse_exp = timeseries(imresize(VSFs1t1_impulse_exp.data,[101 1]),imresize(VSFs1t1_exp.time,[101 1])); %resample time series to 101 data points
VSFs1t2_exp = timeseries(imresize(VSFs1t2_exp.data,[101 1]),imresize(VSFs1t2_exp.time,[101 1])); %resample time series to 101 data points
VSFs1t2_impulse_exp = timeseries(imresize(VSFs1t2_impulse_exp.data,[101 1]),imresize(VSFs1t2_exp.time,[101 1])); %resample time series to 101 data points
VSFs1t3_exp = timeseries(imresize(VSFs1t3_exp.data,[101 1]),imresize(VSFs1t3_exp.time,[101 1])); %resample time series to 101 data points
VSFs1t3_impulse_exp = timeseries(imresize(VSFs1t3_impulse_exp.data,[101 1]),imresize(VSFs1t3_exp.time,[101 1])); %resample time series to 101 data points

VSFs2t1_exp = timeseries(imresize(VSFs2t1_exp.data,[101 1]),imresize(VSFs2t1_exp.time,[101 1])); %resample time series to 101 data points
VSFs2t1_impulse_exp = timeseries(imresize(VSFs2t1_impulse_exp.data,[101 1]),imresize(VSFs2t1_exp.time,[101 1])); %resample time series to 101 data points
VSFs2t2_exp = timeseries(imresize(VSFs2t2_exp.data,[101 1]),imresize(VSFs2t2_exp.time,[101 1])); %resample time series to 101 data points
VSFs2t2_impulse_exp = timeseries(imresize(VSFs2t2_impulse_exp.data,[101 1]),imresize(VSFs2t2_exp.time,[101 1])); %resample time series to 101 data points
VSFs2t3_exp = timeseries(imresize(VSFs2t3_exp.data,[101 1]),imresize(VSFs2t3_exp.time,[101 1])); %resample time series to 101 data points
VSFs2t3_impulse_exp = timeseries(imresize(VSFs2t3_impulse_exp.data,[101 1]),imresize(VSFs2t3_exp.time,[101 1])); %resample time series to 101 data points

VSFs3t1_exp = timeseries(imresize(VSFs3t1_exp.data,[101 1]),imresize(VSFs3t1_exp.time,[101 1])); %resample time series to 101 data points
VSFs3t1_impulse_exp = timeseries(imresize(VSFs3t1_impulse_exp.data,[101 1]),imresize(VSFs3t1_exp.time,[101 1])); %resample time series to 101 data points
VSFs3t2_exp = timeseries(imresize(VSFs3t2_exp.data,[101 1]),imresize(VSFs3t2_exp.time,[101 1])); %resample time series to 101 data points
VSFs3t2_impulse_exp = timeseries(imresize(VSFs3t2_impulse_exp.data,[101 1]),imresize(VSFs3t2_exp.time,[101 1])); %resample time series to 101 data points
VSFs3t3_exp = timeseries(imresize(VSFs3t3_exp.data,[101 1]),imresize(VSFs3t3_exp.time,[101 1])); %resample time series to 101 data points
VSFs3t3_impulse_exp = timeseries(imresize(VSFs3t3_impulse_exp.data,[101 1]),imresize(VSFs3t3_exp.time,[101 1])); %resample time series to 101 data points

%COP

%trial 2
% VSFGRFR_exp_t2 = timeseries((sqrt(t1_exp.forces.FP2(:,3).^2 + t1_exp.forces.FP2(:,2).^2 + ...
%     t1_exp.forces.FP2(:,1).^2)),time); 
% idx1 = find(VSFGRFR_exp_t2.data > 1,1,'first');
% idx2 = find(VSFGRFR_exp_t2.data > 1,1,'last');
% VSFt2_exp = timeseries(VSFGRFR_exp_t2.data(idx1-buffer:idx2+buffer),...
%     VSFGRFR_exp_t2.time(idx1-buffer:idx2+buffer));

%trial 3
% VSFGRFR_exp_t3 = timeseries((sqrt(t1_exp.forces.FP2(:,3).^2 + t1_exp.forces.FP2(:,2).^2 + ...
%     t1_exp.forces.FP2(:,1).^2)),time); 
% idx1 = find(VSFGRFR_exp_t3.data > 1,1,'first');
% idx2 = find(VSFGRFR_exp_t3.data > 1,1,'last');
% VSFt3_exp = timeseries(VSFGRFR_exp_t3.data(idx1-buffer:idx2+buffer),...
%     VSFGRFR_exp_t3.time(idx1-buffer:idx2+buffer));

%% Refit Simulation data to 101 data points
%trial 1  
% VSFGRFR_sim_t1 = timeseries((sqrt(t1_sim.forces.FP2(:,3).^2 + t1_sim.forces.FP2(:,2).^2 + ...
%     t1_sim.forces.FP2(:,1).^2)),time); 
% idx1 = find(VSFGRFR_sim_t1.data > 1,1,'first');
% idx2 = find(VSFGRFR_sim_t1.data > 1,1,'last');
% VSFt1_sim = timeseries(VSFGRFR_sim_t1.data(idx1-buffer:idx2+buffer),...
%     VSFGRFR_sim_t1.time(idx1-buffer:idx2+buffer));
VSFs1t1_sim = timeseries(imresize(VSFs1t1_sim.data,[101 1]),imresize(VSFs1t1_sim.time,[101 1])); %resample time series to 101 data points
VSFs1t1_impulse_sim = timeseries(imresize(VSFs1t1_impulse_sim.data,[101 1]),imresize(VSFs1t1_sim.time,[101 1])); %resample time series to 101 data points
VSFs1t2_sim = timeseries(imresize(VSFs1t2_sim.data,[101 1]),imresize(VSFs1t2_sim.time,[101 1])); %resample time series to 101 data points
VSFs1t2_impulse_sim = timeseries(imresize(VSFs1t2_impulse_sim.data,[101 1]),imresize(VSFs1t2_sim.time,[101 1])); %resample time series to 101 data points
VSFs1t3_sim = timeseries(imresize(VSFs1t3_sim.data,[101 1]),imresize(VSFs1t3_sim.time,[101 1])); %resample time series to 101 data points
VSFs1t3_impulse_sim = timeseries(imresize(VSFs1t3_impulse_sim.data,[101 1]),imresize(VSFs1t3_sim.time,[101 1])); %resample time series to 101 data points

VSFs2t1_sim = timeseries(imresize(VSFs2t1_sim.data,[101 1]),imresize(VSFs2t1_sim.time,[101 1])); %resample time series to 101 data points
VSFs2t1_impulse_sim = timeseries(imresize(VSFs2t1_impulse_sim.data,[101 1]),imresize(VSFs2t1_sim.time,[101 1])); %resample time series to 101 data points
VSFs2t2_sim = timeseries(imresize(VSFs2t2_sim.data,[101 1]),imresize(VSFs2t2_sim.time,[101 1])); %resample time series to 101 data points
VSFs2t2_impulse_sim = timeseries(imresize(VSFs2t2_impulse_sim.data,[101 1]),imresize(VSFs2t2_sim.time,[101 1])); %resample time series to 101 data points
VSFs2t3_sim = timeseries(imresize(VSFs2t3_sim.data,[101 1]),imresize(VSFs2t3_sim.time,[101 1])); %resample time series to 101 data points
VSFs2t3_impulse_sim = timeseries(imresize(VSFs2t3_impulse_sim.data,[101 1]),imresize(VSFs2t3_sim.time,[101 1])); %resample time series to 101 data points

VSFs3t1_sim = timeseries(imresize(VSFs3t1_sim.data,[101 1]),imresize(VSFs3t1_sim.time,[101 1])); %resample time series to 101 data points
VSFs3t1_impulse_sim = timeseries(imresize(VSFs3t1_impulse_sim.data,[101 1]),imresize(VSFs3t1_sim.time,[101 1])); %resample time series to 101 data points
VSFs3t2_sim = timeseries(imresize(VSFs3t2_sim.data,[101 1]),imresize(VSFs3t2_sim.time,[101 1])); %resample time series to 101 data points
VSFs3t2_impulse_sim = timeseries(imresize(VSFs3t2_impulse_sim.data,[101 1]),imresize(VSFs3t2_sim.time,[101 1])); %resample time series to 101 data points
VSFs3t3_sim = timeseries(imresize(VSFs3t3_sim.data,[101 1]),imresize(VSFs3t3_sim.time,[101 1])); %resample time series to 101 data points
VSFs3t3_impulse_sim = timeseries(imresize(VSFs3t3_impulse_sim.data,[101 1]),imresize(VSFs3t3_sim.time,[101 1])); %resample time series to 101 data points


%trial 2
% VSFGRFR_sim_t2 = timeseries((sqrt(t1_sim.forces.FP2(:,3).^2 + t1_sim.forces.FP2(:,2).^2 + ...
%     t1_sim.forces.FP2(:,1).^2)),time); 
% idx1 = find(VSFGRFR_sim_t2.data > 1,1,'first');
% idx2 = find(VSFGRFR_sim_t2.data > 1,1,'last');
% VSFt2_sim = timeseries(VSFGRFR_sim_t2.data(idx1-buffer:idx2+buffer),...
%     VSFGRFR_sim_t2.time(idx1-buffer:idx2+buffer));

%trial 3
% VSFGRFR_sim_t3 = timeseries((sqrt(t1_sim.forces.FP2(:,3).^2 + t1_sim.forces.FP2(:,2).^2 + ...
%     t1_sim.forces.FP2(:,1).^2)),time); 
% idx1 = find(VSFGRFR_sim_t3.data > 1,1,'first');
% idx2 = find(VSFGRFR_sim_t3.data > 1,1,'last');
% VSFt3_sim = timeseries(VSFGRFR_sim_t3.data(idx1-buffer:idx2+buffer),...
%     VSFGRFR_sim_t3.time(idx1-buffer:idx2+buffer));

%%

VSFGRFRs1_exp = [VSFs1t3_exp.data VSFs1t2_exp.data VSFs1t1_exp.data];
VSFGRFRs1_sim = [VSFs1t3_sim.data VSFs1t2_sim.data VSFs1t1_sim.data];
VSFGRFRs1_impulse_exp = [VSFs1t3_impulse_exp.data VSFs1t2_impulse_exp.data VSFs1t1_impulse_exp.data];
VSFGRFRs1_impulse_sim = [VSFs1t3_impulse_sim.data VSFs1t2_impulse_sim.data VSFs1t1_impulse_sim.data];

VSFGRFRs2_exp = [VSFs2t3_exp.data VSFs2t2_exp.data VSFs2t1_exp.data];
VSFGRFRs2_sim = [VSFs2t3_sim.data VSFs2t2_sim.data VSFs2t1_sim.data];
VSFGRFRs2_impulse_exp = [VSFs2t3_impulse_exp.data VSFs2t2_impulse_exp.data VSFs2t1_impulse_exp.data];
VSFGRFRs2_impulse_sim = [VSFs2t3_impulse_sim.data VSFs2t2_impulse_sim.data VSFs2t1_impulse_sim.data];

VSFGRFRs3_exp = [VSFs3t3_exp.data VSFs3t2_exp.data VSFs3t1_exp.data];
VSFGRFRs3_sim = [VSFs3t3_sim.data VSFs3t2_sim.data VSFs3t1_sim.data];
VSFGRFRs3_impulse_exp = [VSFs3t3_impulse_exp.data VSFs3t2_impulse_exp.data VSFs3t1_impulse_exp.data];
VSFGRFRs3_impulse_sim = [VSFs3t3_impulse_sim.data VSFs3t2_impulse_sim.data VSFs3t1_impulse_sim.data];


%Cop
VSFCOPs1_exp = [VSFs1t1_COP_exp.data VSFs1t2_COP_exp.data VSFs1t3_COP_exp.data];
VSFCOPs1_sim = [VSFs1t1_COP_sim.data VSFs1t2_COP_sim.data VSFs1t3_COP_sim.data];
VSFCOPs2_exp = [VSFs2t1_COP_exp.data VSFs2t2_COP_exp.data VSFs2t3_COP_exp.data];
VSFCOPs2_sim = [VSFs2t1_COP_sim.data VSFs2t2_COP_sim.data VSFs2t3_COP_sim.data];
VSFCOPs3_exp = [VSFs3t1_COP_exp.data VSFs3t2_COP_exp.data VSFs3t3_COP_exp.data];
VSFCOPs3_sim = [VSFs3t1_COP_sim.data VSFs3t2_COP_sim.data VSFs3t3_COP_sim.data];

% VSFGRFRS1_exp_mean = mean(VSFGRFRS1_exp,2);
% VSFGRFRS1_exp_sd = std(VSFGRFRS1_exp,[],2);
% SD_upper = (VSFGRFRS1_exp_mean+VSFGRFRS1_exp_sd);
% SD_lower = (VSFGRFRS1_exp_mean-VSFGRFRS1_exp_sd);


subplot(2,3,1);
stdshade(VSFGRFRs1_exp',0.5,'b');hold on
stdshade(VSFGRFRs1_sim',0.5,'r');
axis square;
legend('Experimental SD', 'Experimental mean', 'Simulation SD', 'Simulation mean');
xlabel('Stance phase (percent)','FontWeight','bold','FontSize',18,'FontName','Times');xticks([19 54 95])
xticklabels({'0' '50' '100'})
ylabel('GRF_R (N)','FontWeight','bold','FontSize',18,'FontName','Times');
title('Low Stiffness','FontWeight','Bold','FontSize',24,'FontName','Times'); hold on;
subplot(2,3,4);
% h = figure;
stdshade(VSFGRFRs1_impulse_exp',0.5,'b');hold on
stdshade(VSFGRFRs1_impulse_sim',0.5,'r');hold on
axis square;
legend('Experimental SD', 'Experimental mean', 'Simulation SD', 'Simulation mean');
xlabel('Stance phase (percent)','FontWeight','bold','FontSize',18,'FontName','Times');xticks([19 54 95])
xticklabels({'0' '50' '100'})
ylabel('GRF_R Impulse (N\cdots)','FontWeight','bold','FontSize',18,'FontName','Times');
hold on;
subplot(2,3,2);
stdshade(VSFGRFRs2_exp',0.5,'b');hold on
stdshade(VSFGRFRs2_sim',0.5,'r');
axis square;
legend('Experimental SD', 'Experimental mean', 'Simulation SD', 'Simulation mean');
xlabel('Stance phase (percent)','FontWeight','bold','FontSize',18,'FontName','Times');xticks([19 54 95])
xticklabels({'0' '50' '100'})
ylabel('GRF_R (N)','FontWeight','bold','FontSize',18,'FontName','Times');
title('Medium Stiffness','FontWeight','Bold','FontSize',24,'FontName','Times'); hold on;
subplot(2,3,5);
% h = figure;
stdshade(VSFGRFRs2_impulse_exp',0.5,'b');hold on
stdshade(VSFGRFRs2_impulse_sim',0.5,'r');hold on
axis square;
legend('Experimental SD', 'Experimental mean', 'Simulation SD', 'Simulation mean');
xlabel('Stance phase (percent)','FontWeight','bold','FontSize',18,'FontName','Times');xticks([19 54 95])
xticklabels({'0' '50' '100'})
ylabel('GRF_R Impulse (N\cdots)','FontWeight','bold','FontSize',18,'FontName','Times');
hold on;
hold on;
subplot(2,3,3);
stdshade(VSFGRFRs3_exp',0.5,'b');hold on
stdshade(VSFGRFRs3_sim',0.5,'r');
axis square;
legend('Experimental SD', 'Experimental mean', 'Simulation SD', 'Simulation mean');
xlabel('Stance phase (percent)','FontWeight','bold','FontSize',18,'FontName','Times');xticks([19 54 95])
xticklabels({'0' '50' '100'})
ylabel('GRF_R (N)','FontWeight','bold','FontSize',18,'FontName','Times');
title('High Stiffness','FontWeight','Bold','FontSize',24,'FontName','Times'); hold on;
subplot(2,3,6);
% h = figure;
stdshade(VSFGRFRs3_impulse_exp',0.5,'b');hold on
stdshade(VSFGRFRs3_impulse_sim',0.5,'r');hold on
axis square;
legend('Experimental SD', 'Experimental mean', 'Simulation SD', 'Simulation mean');
xlabel('Stance phase (percent)','FontWeight','bold','FontSize',18,'FontName','Times');xticks([19 54 95])
xticklabels({'0' '50' '100'})
ylabel('GRF_R Impulse (N\cdots)','FontWeight','bold','FontSize',18,'FontName','Times');
sgtitle('Simulation vs Experimental GRFR (top) and Impulse (Bottom)')
hold off;

%%
figure; hold on;
subplot(1,3,1);
stdshade(VSFCOPs1_exp',0.5,'b');hold on
stdshade(VSFCOPs1_sim',0.5,'r');
axis square;
legend('Experimental SD', 'Experimental mean', 'Simulation SD', 'Simulation mean');
xlabel('Stance phase (percent)','FontWeight','bold','FontSize',18,'FontName','Times');xticks([0 50.5 101])
xticklabels({'0' '50' '100'})
ylabel('Anterior-posterior COP (percent foot length)','FontWeight','bold','FontSize',18,'FontName','Times');yticks([0 57.25 114.5 171.75 229])
yticklabels({'0' '25' '50' '75' '100'});
title('Low Stiffness','FontWeight','Bold','FontSize',24,'FontName','Times'); hold on;

subplot(1,3,2);
stdshade(VSFCOPs2_exp',0.5,'b');hold on
stdshade(VSFCOPs2_sim',0.5,'r');
axis square;
legend('Experimental SD', 'Experimental mean', 'Simulation SD', 'Simulation mean');
xlabel('Stance phase (percent)','FontWeight','bold','FontSize',18,'FontName','Times');xticks([0 50.5 101])
xticklabels({'0' '50' '100'})
ylabel('Anterior-posterior COP (percent foot length)','FontWeight','bold','FontSize',18,'FontName','Times');yticks([0 57.25 114.5 171.75 229])
yticklabels({'0' '25' '50' '75' '100'});
title('Medium Stiffness','FontWeight','Bold','FontSize',24,'FontName','Times'); hold on;

subplot(1,3,3);
stdshade(VSFCOPs3_exp',0.5,'b');hold on
stdshade(VSFCOPs3_sim',0.5,'r');
axis square;
legend('Experimental SD', 'Experimental mean', 'Simulation SD', 'Simulation mean');
xlabel('Stance phase (percent)','FontWeight','bold','FontSize',18,'FontName','Times');xticks([0 50.5 101])
xticklabels({'0' '50' '100'})
ylabel('Anterior-posterior COP (percent foot length)','FontWeight','bold','FontSize',18,'FontName','Times');yticks([0 57.25 114.5 171.75 229])
yticklabels({'0' '25' '50' '75' '100'});
title('High Stiffness','FontWeight','Bold','FontSize',24,'FontName','Times'); 
