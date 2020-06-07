%% Parameterize model
%Gravity
if isequal(GravityOn,1) 
    gravity = -9.80665;
else
    gravity = 0;
end

%% Walker Model GEOMETRY

%variable stiffness foot geometry
keel_length_total = 229; %mm
keel_thickness = 6.35; %mm
keel_width = 68.4; %mm
keel_length_rigid = 66; %mm
keel_length_tip = 78; %mm
keel_discretized_length = 10.5; %(keel_length_total-(keel_length_tip+keel_length_rigid))/keel_descretization; %mm
keel_density = 6572.37; %kg/m^3

%General parameters
bone_color = [1 1 1];
density = 1000;
foot_density = 2000;
world_damping = 300; %0.25;
world_rot_damping = 1000; %0.25;
if ~exist('actuatorType','var')
    actuatorType = 1;
end

%ANATOMIC PARAMETERS

%Segment Length parameters
%Foot
foot_x = subj1.static.segLength.footLX ;
foot_y = subj1.static.segLength.footLY;
foot_z = subj1.static.segLength.footLZ;
    %Foot extrusion
        % w = foot_z;
        % l = foot_x-2;
        % A = linspace(-pi/2, pi/2, 20)';
        % B = linspace(pi/2, 3*pi/2, 20)';
        % csRight = [l/4 + w/1.2*cos(A) w/4*sin(A)];
        % csLeft = [-l/4 + w/1.2*cos(B) w/4*sin(B)];
        % cs = [csRight; csLeft];
        % FootExtrusion_opc = 0.2;
        % Plot foot extrusion
        % figure; hold on; axis equal;
        % plot(cs(:,1), cs(:,2), 'Color', [0.6 0.6 0.6], 'Marker', '.',...
        % 'MarkerSize', 9, 'MarkerEdgeColor', [1 0 0]); 

osf = 1; %Offset scale factor (used to adjust joint and segment offsets when scaling between 

%Leg
leg_radius = 9.5;
shank_radius = 5.5;

lower_leg_length = subj1.static.segLength.shankLengthL;
upper_leg_length = subj1.static.segLength.legLengthL;
lower_leg_length_amputated = subj1.static.segLength.shankLengthR;
res_shank_length = lower_leg_length_amputated/lower_leg_length; %residual leg length expressed as a % of intact shank length
pylon_radius = 1.5; %cm
socket2shank_offset = [1.5 1 25.5]; %cm 
pylon_length = lower_leg_length - lower_leg_length_amputated - (socket2shank_offset(3)/10) + 8.7; %cm
pylon_density = 5000; %kg/m^3

%Torso (cm)
torso_x = subj1.static.segLength.pelvisX;
torso_y = subj1.static.segLength.pelvisY;
torso_z = subj1.static.segLength.pelvisZ;
% torso_offset_z = -.1;
% torso_offset_x = -0.5;
% torso_offset_y = .1;

%Mass parameters
SubjMass = subj1.static.mass; %kg
mass_HAT = 0.628*SubjMass; % regression model limb masses from de Leva, 1996
mass_leg = .1416*SubjMass;
    volume_leg = ((pi/3)*(upper_leg_length/100))*((leg_radius/100)^2+(shank_radius/100)^2+(leg_radius/100*shank_radius/100));
    leg_density = mass_leg/volume_leg;
mass_shank = .0433*SubjMass;
mass_shank_amputated = mass_shank*res_shank_length;
    volume_shank = ((pi/3)*(lower_leg_length/100))*((shank_radius/100)^2+(shank_radius/100)^2+(shank_radius/100*shank_radius/100));
    shank_density = mass_shank/volume_shank;

mass_foot = .0137*SubjMass;
mass_total = mass_HAT+mass_leg*2+mass_shank+mass_shank_amputated+mass_foot+.08;
mass_pylon = .5*(SubjMass-mass_total); %kg
mass_socket = .5*(SubjMass-mass_total); %kg
mass_total = mass_HAT+mass_leg*2+mass_socket+mass_pylon+mass_shank+mass_shank_amputated+mass_foot+.08;

%% Joint offsets %socket to shank offset listed under segment length parameters

RHip2Torso_offset = subj1.static.offsets.LegRProxPos_pelvis; %cm %%TRY GETTING PELVIS POS IN GCS AND SUBTRACTING 
LHip2Torso_offset = subj1.static.offsets.LegLProxPos_pelvis;
    y = [1 0 0;...
        0 1 0;...
        0 0 -1];
    z = [1 0 0;...
        0 1 0;...
        0 0 -1];
RKnee2Leg_offset = y*subj1.static.offsets.FKneeRPos_LegR; clear y
LKnee2Leg_offset = z*subj1.static.offsets.FKneeLPos_LegL; clear z
LShank2Knee_offset = osf*[0 0 2];
RShank2Knee_offset = osf*[0 0 2];
LShank2Ankle_offset = subj1.static.offsets.FAJCLPos_ShankL;

% 
% RHip2Torso_offset = osf*[2.70 -6.40 -9.16]; %cm
% LHip2Torso_offset = osf*[2.70 6.40 -9.16];
% RKnee2Leg_offset = osf*[-.2 subj1.seg.legLengthR -1.10]; %Expressed in knee coordinate system (rotated -90 about X wrt leg coordinate system)
% LKnee2Leg_offset = osf*[-.2 subj1.seg.legLengthL 1.10]; 
% LShank2Knee_offset = osf*[0 0 2];
% RShank2Knee_offset = osf*[0 0 2];
% LShank2Ankle_offset = osf*[2 0 subj1.seg.shankLengthL+.75];
foot_offset_vsf = [2.5 1.5 -3];
foot_offset = [-2.5 -1.5 4.5];
footL_anatomic_offset = [-1 .5 -.8];
% % RHip2Torso_offset = [.09*subj1.seg.ASISYdist -.33*subj1.seg.ASISYdist -.35*subj1.seg.ASISYdist]; %cm
% LHip2Torso_offset = [.09*subj1.seg.ASISYdist .33*subj1.seg.ASISYdist -.35*subj1.seg.ASISYdist]; 
% RKnee2Leg_offset = [subj1.seg.legLengthL*.05 subj1.seg.legLengthL+2 subj1.seg.legLengthL*.02]; %Change to Right to account for assymentries (future adaptation)
% LKnee2Leg_offset = [subj1.seg.legLengthL*.05 subj1.seg.legLengthL+1.5 -subj1.seg.legLengthL*.01];
% LShank2Ankle_offset = [0 .025*subj1.seg.shankLengthL subj1.seg.shankLengthL];
% 


%Contact parameters
contact_stiffness = 5000;
contact_damping =.7*contact_stiffness;
mu_k = 0.99;
mu_s = 0.99;
mu_vth = 0.01;
height_plane = 0.025;
plane_x = 25;
plane_y = 3;
contact_point_radius = .001;
contact_color = [0 0 0];
contact_opc = 0;

init_height = foot_z + subj1.static.segLength.legLengthL + subj1.static.segLength.shankLengthL + ...
    subj1.static.segLength.pelvisZ + height_plane/2;

%%

%Contact/friction parameters DEFAULTS
% contact_stiffness = 2500;
% contact_damping = 100;
% mu_k = 0.6;
% mu_s = 0.8;
% mu_vth = 0.1;
% height_plane = 0.025; %m
% plane_x = 25;
% plane_y = 3;
% contact_point_radius = .01;
% contact_color = [0 0 0];
% contact_opc = 0;
% 
% %Foot parameters (cm)
% foot_x = 22;
% foot_y = 7;
% foot_z = 1;
% foot_offset_vsf = [-2.7 31 -3.3];
% % foot_offset = [-1 0 0];
% foot_offset = [-7 0 0];
% 
% %Leg parameters (cm)
% 
% leg_radius = 2;
% lower_leg_length = 45.5;
% upper_leg_length = 40.5;
% lower_leg_length_amputated = 10;
% 
% %Torso parameters (cm)
% torso_y = 20;
% torso_x = 10;
% torso_z = 16;
% torso_offset_z = -2;
% torso_offset_x = -0.5;
% init_height = foot_z + lower_leg_length + upper_leg_length + ...
%               torso_z/2 + torso_offset_z + height_plane/2;

%Joint parameters
joint_damping = 0;
joint_stiffness = 0;
motion_time_constant = (60/subj1.framerate{1})/2; %0.025;

%Joint controller parameters
hip_servo_kp = 60;
hip_servo_ki = 10;
hip_servo_kd = 20;
knee_servo_kp = 60;
knee_servo_ki = 5;
knee_servo_kd = 10;
ankle_servo_kp = 20;
ankle_servo_ki = 4;
ankle_servo_kd = 8;
deriv_filter_coeff = 100;
max_torque = 20;

%electric motor parameters
motor_resistance = 1;
motor_constant = 0.02;
motor_inertia = 0;
motor_damping = 0;
motor_inductance = 1.2e-6;
gear_ratio = 50;


%% VSF DYNAMICS
%Revolute joint parameters
revolute_equilibrium_position = 0; %deg
revolute_position_target = 0; %deg
k_inf = 1000000;%N*m/deg
krot = (-0.0677*fulcrum_position)+13.839; %Regression equation for FP vs k (R^2 > .997)

%set up rotational stiffness variables
revolute_66_stiffness = krot;
revolute_76_5_stiffness = k_inf;
revolute_87_stiffness = krot;
revolute_97_5_stiffness = k_inf;
revolute_108_stiffness = krot;
revolute_118_5_stiffness = k_inf;
revolute_129_stiffness = krot;
revolute_139_5_stiffness = k_inf;
revolute_150_stiffness = krot; 
revolute_160_5_stiffness = k_inf; 
revolute_171_stiffness = krot;
revolute_181_5_stiffness = k_inf;
revolute_192_stiffness = krot;
revolute_202_5_stiffness = k_inf;
revolute_213_stiffness = krot;

if isequal(fulcrum_position,87)
    revolute_66_stiffness = k_inf; 
    revolute_76_5_stiffness = k_inf;
    revolute_87_stiffness = 6.75;
    revolute_97_5_stiffness = k_inf;
    revolute_108_stiffness = 6.75;
    revolute_118_5_stiffness = k_inf;
    revolute_129_stiffness = 6.75;
    revolute_139_5_stiffness = k_inf;
    revolute_150_stiffness = 6.75; 
    revolute_160_5_stiffness = k_inf; 
    revolute_171_stiffness = 6.75;
    revolute_181_5_stiffness = k_inf;
    revolute_192_stiffness = 6.75;
    revolute_202_5_stiffness = k_inf;
    revolute_213_stiffness = 6.75;
end

if isequal(fulcrum_position,108)
    revolute_66_stiffness = k_inf; 
    revolute_76_5_stiffness = k_inf;
    revolute_87_stiffness = k_inf;
    revolute_97_5_stiffness = k_inf;
    revolute_108_stiffness = krot;
    revolute_118_5_stiffness = k_inf;
    revolute_129_stiffness = krot;
    revolute_139_5_stiffness = k_inf;
    revolute_150_stiffness = krot; 
    revolute_160_5_stiffness = k_inf; 
    revolute_171_stiffness = krot;
    revolute_181_5_stiffness = k_inf;
    revolute_192_stiffness = krot;
    revolute_202_5_stiffness = k_inf;
    revolute_213_stiffness = krot;
end

if isequal(fulcrum_position,129)
    revolute_66_stiffness = k_inf; 
    revolute_76_5_stiffness = k_inf;
    revolute_87_stiffness = k_inf;
    revolute_97_5_stiffness = k_inf;
    revolute_108_stiffness = k_inf;
    revolute_118_5_stiffness = k_inf;
    revolute_129_stiffness = krot;
    revolute_139_5_stiffness = k_inf;
    revolute_150_stiffness = krot; 
    revolute_160_5_stiffness = k_inf; 
    revolute_171_stiffness = krot;
    revolute_181_5_stiffness = k_inf;
    revolute_192_stiffness = krot;
    revolute_202_5_stiffness = k_inf;
    revolute_213_stiffness = krot;
end

if isequal(fulcrum_position,150)
    revolute_66_stiffness = k_inf; 
    revolute_76_5_stiffness = k_inf;
    revolute_87_stiffness = k_inf;
    revolute_97_5_stiffness = k_inf;
    revolute_108_stiffness = k_inf;
    revolute_118_5_stiffness = k_inf;
    revolute_129_stiffness = k_inf;
    revolute_139_5_stiffness = k_inf;
    revolute_150_stiffness = krot; 
    revolute_160_5_stiffness = k_inf; 
    revolute_171_stiffness = krot;
    revolute_181_5_stiffness = k_inf;
    revolute_192_stiffness = krot;
    revolute_202_5_stiffness = k_inf;
    revolute_213_stiffness = krot;
end

%Set damping factor (N*m(deg/s))
if isequal(fulcrum_position,150) 
   revolute_damping = krot*.001;
else 
   revolute_damping = krot*.5;
end

%% Motion Inputs
if isequal(ModelControl,0)
    gaitPeriod = 0.8;
    time = linspace(0,gaitPeriod,7)';
    ankle_motion = deg2rad([-7.5 10 10 5 0 -10 -7.5]');
    knee_motion = deg2rad([10, -5, 2.5, -10, -10, 15, 10]');
    hip_motion = deg2rad([-10, -7.5, -15, 10, 15, 10, -10]');
    curveDataL = createSmoothTrajectory(ankle_motion,knee_motion,hip_motion,gaitPeriod);
    curveDataR = createSmoothTrajectory(ankle_motion,knee_motion,hip_motion,gaitPeriod);
    timeDelay = .5*gaitPeriod;
    simTime = 10;
else 
    gaitPeriod = length(subj1.km.Rankle(:,1))/subj1.framerate{1};
    time = linspace(0,gaitPeriod,length(subj1.km.Rankle))';
    
    pelvisPosX = [time,subj1.km.PelvisCOMPos(:,1)];
    pelvisPosY = [time,subj1.km.PelvisCOMPos(:,2)];  
    pelvisPosZ = [time,subj1.km.PelvisCOMPos(:,3)];
    pelvisAngX = [time,deg2rad(subj1.km.PelvisAngle(:,1))];
    pelvisAngY = [time,deg2rad(subj1.km.PelvisAngle(:,2))];
    pelvisAngZ = [time,deg2rad(subj1.km.PelvisAngle(:,3))];
    
    hipRAngX = [time,deg2rad(subj1.km.Rhip(:,1))];
    hipRAngY = [time,-1*deg2rad(subj1.km.Rhip(:,2))];
    hipRAngZ = [time,deg2rad(subj1.km.Rhip(:,3))]; 
    hipLAngX = [time,-1*deg2rad(subj1.km.Lhip(:,1))];
    hipLAngY = [time,-1*deg2rad(subj1.km.Lhip(:,2))];
    hipLAngZ = [time,-1*deg2rad(subj1.km.Lhip(:,3))]; 
    
    ankle_motionRX = [time,deg2rad(subj1.km.Rankle(:,1))];
        a = linspace(0,0,length(time));
    ankle_motionRY = [time,a']; clear a
    ankle_motionRZ = [time,deg2rad(subj1.km.Rankle(:,3))];
    ankle_motionLX = [time,-1*deg2rad(subj1.km.Lankle(:,1))];
    ankle_motionLY = [time,-1*deg2rad(subj1.km.Lankle(:,2))];
    ankle_motionLZ = [time,-1*deg2rad(subj1.km.Lankle(:,3))];

    knee_motionL = [time,deg2rad(subj1.km.Lknee(:,2))];
    knee_motionR = [time,deg2rad(subj1.km.Rknee(:,2))];
    
    timeDelay = 0;
    simTime = gaitPeriod;
end


%% 

%% Contact model parameters
contact_stiffness = 40000;
contact_damping =.6*contact_stiffness;

contact_stiffness_heel = 8000;
contact_damping_heel =.85*contact_stiffness_heel;

contact_stiffness_intact = 15000;
contact_damping_intact =.6*contact_stiffness_intact;


mu_k = .999;
mu_s = 0.99;
mu_vth = .1;
height_plane = 0.025;
plane_x = 25;
plane_y = 3;
contact_point_radius = .001;
contact_color = [1 0 0];
contact_opc = .0;


contact_color_intact = [1 1 0];
contact_opc_intact = 0;

init_height = foot_z + subj1.static.segLength.legLengthL + subj1.static.segLength.shankLengthL + ...
    subj1.static.segLength.pelvisZ + height_plane/2;


world_model_damping = 1;
world_model_stiffness = 1;

bone_color = [1 1 1];

limb_geom_opc = 0;

%% Open model
% GaitModel_LPVSF_R_BodyParts3D

GaitModel_LPVSF_R_wscv3d_simpintactCM_2segfoot_adjustLCS;

%% Run simulations and compare the results with experimental data
% sim('GaitModel_LPVSF_R_BodyParts3D');
% 
% CompareSim

%% MTS parameters
%MTS geometry
% mts_OuterRadius = .01; %m
% mts_length = .068; %m
% mts_density = .001; %kg/m^3
% 
% %MTS kinematics
% %for motion control (comment out for force control):
% max_disp = 26.894; %mm
% mts_displacement = linspace(42.224,42.224+max_disp,100); %[5;20;35;50;65;80;95]; %mm
% mts_contact_damping = 100; %N/m/s
% mts_contact_stiffness = 1000000; %N/m
% %for force control:
% mts_force = linspace(0,800,100);

% %% Bell and brand HJC transforms
% RPelvis2HJC_offset = [10.24, -14.5, -1.95] + [-.3*29, .14*29, -.22*29];
% LPelvis2HJC_offset = [10.24, 14.5, -1.95] + [-.3*29, -.14*29, -.22*29];

% HJCL = [LASIavg - (.3*ASISdist), LASIavg + (.14*ASISdist), LASIavg - (.22*ASISdist)]; %Bell and Brand regression model
