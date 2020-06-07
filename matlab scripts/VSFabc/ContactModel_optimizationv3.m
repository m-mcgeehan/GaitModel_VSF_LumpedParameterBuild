
FP2_trim = dynamic.forces.FP2(200:475,3); %VSF side
FP2time_trim = time(200:475);

FP2_DS = downsample(FP2_trim,2);
FP2time_DS = downsample(FP2time_trim,2);


FP1_trim = dynamic.forces.FP1(160:350,3); %Intact side
FP1time_trim = time(160:350); 

FP1_DS = downsample(FP1_trim,2);
FP1time_DS = downsample(FP1time_trim,2);

%%
FP2Resultant_trim = sqrt(dynamic.forces.FP2(200:475,3).^2 + dynamic.forces.FP2(200:475,2).^2 + ...
    dynamic.forces.FP2(200:475,1).^2); %VSF side
FP2time_trim = time(200:475);

FP2Resultant_DS = downsample(FP2Resultant_trim,2);
FP2time_DS = downsample(FP2time_trim,2);


FP1Resultant_trim = sqrt(dynamic.forces.FP1(160:350,3).^2 + dynamic.forces.FP1(160:350,2).^2 ...
    + dynamic.forces.FP1(160:350,1).^2); %Intact side
FP1time_trim = time(160:350); 

FP1Resultant_DS = downsample(FP1Resultant_trim,2);
FP1time_DS = downsample(FP1time_trim,2);

subplot(1,2,1); plot(FP2time_DS, FP2Resultant_DS,'k'); hold on;
plot(FP2time_DS, FP2_DS,'r'); legend('Resultant','Vertical'); title('VSF GRF')
subplot(1,2,2); plot(FP1time_DS, FP1Resultant_DS,'k'); hold on;
plot(FP1time_DS, FP1_DS,'r'); legend('Resultant','Vertical'); title('Intact GRF')


set_param('GaitModelVSF_2_24_20','StartTime','0','StopTime','2.3689'); %Intact and VSF GRF

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

heel_prisK = 18000;
heel_prisD = 10000;
centerX = 6;
centerZ = 2;
%% For HCM2 heel contact model
ck_heel_VSF = 50000;
cd_heel_VSF = 500;
ck_heel_VSF_post = 50000;
cd_heel_VSF_post = 500;
penDepth_heel_VSF = .01;
penExp_heel_VSF = 1.1;
penDepth_heel_VSF_post = .01;
penExp_heel_VSF_post = 1.1;
rad_heel_VSF = .03;
ztrans_heel_VSF = 3;

VSF_heel_prisK = 8000;
VSF_heel_prisD = 7000; 

pyrK = 2000;
pyrD = 1000;


%% Parameterize model from SDO
%Transfer data from simulink design optimization, unpack parm.continuous,
%and overwrite workspace variables

for i = 1:length(VSF_designVars)
    var = VSF_designVars(i,1);
    varName = var.Name;
    varValue = var.Value;
    designvars.(varName) = varValue;
end
clear i var varName varValue

struct2vars(designvars);

cd_heel_VSF = designvars.cd_heel_VSF;                                                 
ck_heel_VSF = designvars.ck_heel_VSF;                                                 
contact_damping_intact_MTP = designvars.contact_damping_intact_MTP;                   
contact_damping_intact_forefoot = designvars.contact_damping_intact_forefoot;         
contact_damping_intact_heel = designvars.contact_damping_intact_heel;                 
contact_damping_intact_lateral = designvars.contact_damping_intact_lateral;           
contact_damping = designvars.contact_damping;                                         
contact_point_radius = designvars.contact_point_radius;                               
contact_stiffness = designvars.contact_stiffness;                                     
contact_stiffness_intact_MTP = designvars.contact_stiffness_intact_MTP;               
contact_stiffness_intact_forefoot = designvars.contact_stiffness_intact_forefoot;     
contact_stiffness_intact_heel = designvars.contact_stiffness_intact_heel;             
contact_stiffness_intact_lateral = designvars.contact_stiffness_intact_lateral;       
heel_intact_rad = designvars.heel_intact_rad;                                         
heel_prisD = designvars.heel_prisD;                                                   
heel_prisK = designvars.heel_prisK;                                                   
penDepth_heel_VSF = designvars.penDepth_heel_VSF;                                     
penExp_heel_VSF = designvars.penExp_heel_VSF;                                         
pen_depth_heel_intact = designvars.pen_depth_heel_intact;                             
pen_depth_lateral_intact = designvars.pen_depth_lateral_intact;                       
pen_exp_heel_intact = designvars.pen_exp_heel_intact;                                 
pen_exp_lateral_intact = designvars.pen_exp_lateral_intact;                           
pyrD = designvars.pyrD;                                                               
pyrK = designvars.pyrK;                                                               
rad_heel_VSF = designvars.rad_heel_VSF;                                               
xtrans = designvars.xtrans;                                                           
ztrans_heel_VSF = designvars.ztrans_heel_VSF;                                         
VSF_heel_prisD = designvars.VSF_heel_prisD;                                           
VSF_heel_prisK = designvars.VSF_heel_prisK;                                           
lateral_intact_rad = designvars.lateral_intact_rad;                                   
socketD_rot = designvars.socketD_rot;                                                 
socketD_trans = designvars.socketD_trans;                                             
socketK_rot = designvars.socketK_rot;                                                 
socketK_trans = designvars.socketK_trans;                                             
mu_k = designvars.mu_k;                                                               
mu_s = designvars.mu_s;                                                               
mu_vth = designvars.mu_vth;                                                           
intact_fc = designvars.intact_fc;                                                     
vsf_fc = designvars.vsf_fc;

%% from optimization 2-21-20
% cd_heel_VSF = 371808;
% ck_heel_VSF = -19589;
% 
% contact_damping_intact_MTP = 78.8;
% contact_damping_intact_forefoot = 204.2;
% contact_damping_intact_heel = 5351.0;
% contact_damping_intact_lateral = 4574.6;
% contact_stiffness_intact_MTP = 4221.0;
% contact_stiffness_intact_forefoot = 13512.8;
% contact_stiffness_intact_heel = 3181.7;
% contact_stiffness_intact_lateral = 6548.6;
% 
% penDepth_heel_VSF = .0014;
% penExp_heel_VSF = 8.07;
% pen_dept_heel_intact = .0225;
% pen_exp_heel_intact = 14.7;
% pen_exp_lateral_intact = 15.8;
% 
% lateral_intact_rad = .025;
% heel_intact_rad = .03;
% 
% heel_prisD = 5930.8;
% heel_prisK = 59519.8;
% 
% centerX = 6;
% centerZ = 2;
% 
% contact_damping = 327.5;
% contact_stiffness = 7153.0;
% 
% mu_k = .8;
% mu_s = .2;
% 
% VSF_heel_prisD = 10286.0;
% VSF_heel_prisK = 17859.0;
% 
% cd_heel_VSF_post = 807.4;
% ck_heel_VSF_post = 89628.9;
% 
% penDepth_heel_VSF_post = .006;
% penExp_heel_VSF_post = .431;
% 
% pyrD = 1312.8;
% pyrK = 2303.7;
% 
% socketK_trans = 20000;
% socketD_trans = 10000;
% socketK_rot = 10000;
% socketD_rot = 5000;
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



