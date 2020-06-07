
FP2_trim = subj1.forces.FP2(200:475,3); %VSF side
FP2time_trim = time(200:475);

FP2_DS = downsample(FP2_trim,2);
FP2time_DS = downsample(FP2time_trim,2);


FP1_trim = subj1.forces.FP1(160:350,3); %Intact side
FP1time_trim = time(160:350); 

FP1_DS = downsample(FP1_trim,2);
FP1time_DS = downsample(FP1time_trim,2);

set_param('HCM_2_GaitModelVSF_2_13_20','StartTime','0','StopTime','2.3689'); %Intact and VSF GRF

% set_param('HCM_2_GaitModelVSF_1_21_20','StartTime','0.996631147540984','StopTime','2.3689'); %VSF GRF

% set_param('HCM_1_GaitModelVSF_1_21_19','StartTime','1.3172','StopTime','2.3689'); %VSF GRF
% set_param('HCM_1_GaitModelVSF_1_21_19','StartTime','.7963','StopTime','1.7479'); %Intact GRF
% set_param('HCM_1_GaitModelVSF_1_21_19','StartTime','0','StopTime','3.0550');

% set_param('Copy_of_GaitModel_VSF_R_11_18_19','StartTime','1.3172','StopTime','2.3689');
% set_param('Copy_of_GaitModel_VSF_R_11_18_19','StartTime','1.3172','StopTime','1.6777');

contact_point_radius_heel = .04;


%% 
contact_d_heel_mid = 1500;
contact_k_heel_mid = 2000;
contact_k_heel_ant = 2000;
contact_d_heel_ant = 1000;
contact_k_heel_post = 5000;
contact_d_heel_post = 4500;
contact_k_heel_pc = 10000;
contact_d_heel_pc = 8000;

pen_exp_heel_intact = 1;
pen_depth_heel_intact = 0.1;
heel_intact_rad = .03;
pen_exp_lateral_intact = 1;
pen_depth_lateral_intact = 0.1;
lateral_intact_rad = .03;
%% from LHC optimzation
contact_k_heel_post = 520;
contact_d_heel_post = 1748;
contact_k_heel_pc = 448;
contact_d_heel_pc = 720;
%% 
contact_r_heel = .001;
xtrans = 12;
ztrans = 1.5;
heel_rad = 0.03;
heel_pen = .01;
pen_exp = 1;

%% From intact optimization
contact_damping_intact_MTP = 118.7;
contact_damping_intact_forefoot = 160.9;
contact_damping_intact_heel = 2910.9;
contact_damping_intact_lateral = 4036.5;

% contact_stiffness_intact_MTP = 118.7;
contact_stiffness_intact_MTP = 10118.7;

% contact_stiffness_intact_forefoot = 160.9;
contact_stiffness_intact_forefoot = 15600.9;

contact_stiffness_intact_heel = 2910.9;
contact_stiffness_intact_lateral = 5874.4;

pen_depth_heel_intact = .00751;
pen_depth_lateral_intact = .00574;
pen_exp_heel_intact = 1.4297;
pen_exp_lateral_intact = 1.3741;
heel_intact_rad = .03;
lateral_intact_rad = .03;

%% For HCM2 heel contact model
ck_heel_VSF = 50000;
cd_heel_VSF = 500;
penDepth_heel_VSF = .01;
penExp_heel_VSF = 1.1;

rad_heel_VSF = .03;
ztrans_heel_VSF = 3;

%% 
simTime = x{3};
forceTime = x{4};
simForce = y{3};
force = y{4};

simForce_interp = interp1(simTime,simForce,forceTime);


order = 2;
fs = 200;
fc = 90;%50 hz cut off
Wn_fp = 2*fc/fs; 
[bfp,afp] = butter(order,Wn_fp); 
simForceR = filtfilt(bfp,afp,simForce_interp(2:end)); %symmetric lowpass filter of normal force



plot(forceTime(2:end),simForceR); hold on;
plot(forceTime,force);


% sampleTime = 60/200;
% endTime = gaitPeriod;
% numberOfSamples = endTime*1/sampleTime +1;
% timeVector = (0:numberOfSamples)*sampleTime;
% 
% pelvisAngX_TS = timeseries(pelvisAngX(:,2),pelvisAngX(:,1));
% 
% t132 = find(time == 1.32);
% t232 = find(time == 1.84);

plot(FP1TS.time,FP1TS.data(:,3)); hold on; plot(FP2TS.time,FP2TS.data(:,3)); hold on; plot(VSF_CF); hold on; plot(intact_CF);



