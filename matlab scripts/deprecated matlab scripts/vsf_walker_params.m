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
density = 1000;
foot_density = 2000;
world_damping = 300; %0.25;
world_rot_damping = 1000; %0.25;
if ~exist('actuatorType','var')
    actuatorType = 1;
end

%ANATOMIC PARAMETERS
%Mass parameters
res_shank_length = .3; %residual leg length expressed as a % of intact shank length
mass_HAT = 0.628*SubjMass; % regression model limb masses from de Leva, 1996
mass_leg = .1416*SubjMass;
mass_pylon = 1.5; %kg
mass_socket = 1.5; %kg
mass_shank = .0433*SubjMass;
mass_shank_amputated = mass_shank*res_shank_length;
mass_foot = .0137*SubjMass;
mass_total = mass_HAT+mass_leg*2+mass_socket+mass_pylon+mass_shank+mass_shank_amputated+mass_foot+.08;

%Length parameters
%Foot
foot_x = 22;
foot_y = 7;
foot_z = 7;
foot_offset_vsf = [2.5 1.5 -3];
% foot_offset = [-1 0 0];
foot_offset = [-7 0 0];

%Leg parameters (cm)
leg_radius = 5.5;
lower_leg_length = .232*SubjHeight;
upper_leg_length = .247*SubjHeight;
lower_leg_length_amputated = lower_leg_length*res_shank_length;
socket_limb_offset = [0 0 8.5]; %cm
pylon_radius = 1.5; %cm
pylon_length = lower_leg_length - lower_leg_length_amputated - socket_limb_offset(1)-9.124; %cm
pylon_density = 50000; %kg/m^3

%Torso parameters (cm)
torso_y = 27;
torso_x = 16;
torso_z = 16;
torso_offset_z = -2;
torso_offset_x = -0.5;


% contact_stiffness = 100000;
% contact_damping = 5000;
contact_stiffness = 5000;
contact_damping =.7*contact_stiffness;
mu_k = 0.99;
mu_s = 0.99;
mu_vth = 0.01;
height_plane = 0.025;
plane_x = 25;
plane_y = 3;
contact_point_radius = .01;
contact_color = [0 0 0];
contact_opc = 0;

init_height = foot_z + lower_leg_length + upper_leg_length + ...
              torso_z/2 + torso_offset_z + height_plane/2;
init_height = init_height - 3.5;
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
motion_time_constant = 0.01; %0.025;

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

%% from robot model


% Motion Inputs
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
    ankle_motionR = deg2rad(subj1.km.Rankle(:,2));
    knee_motionR = deg2rad(subj1.km.Rknee(:,2));
    hip_motionR = deg2rad(subj1.km.Rhip(:,2));
    curveDataR = createSmoothTrajectory(ankle_motionR,knee_motionR,hip_motionR,gaitPeriod);
    ankle_motionL = deg2rad(subj1.km.Lankle(:,2));
    knee_motionL = deg2rad(subj1.km.Lknee(:,2));
    hip_motionL = deg2rad(subj1.km.Lhip(:,2));
    curveDataL = createSmoothTrajectory(ankle_motionL,knee_motionL,hip_motionL,gaitPeriod);
    pelvisPosX = subj1.km.PelvisCOMPos(:,1);
    pelvisPosY = subj1.km.PelvisCOMPos(:,2);
    pelvisPosZ = subj1.km.PelvisCOMPos(:,3);
    pelvisAngX = deg2rad(subj1.km.PelvisAngle(:,1));
    pelvisAngY = deg2rad(subj1.km.PelvisAngle(:,2));
    pelvisAngZ = deg2rad(subj1.km.PelvisAngle(:,3));
    timeDelay = 0;
    simTime = gaitPeriod;
end

pelvisPosY = [time,pelvisPosY];
pelvisPosX = [time,pelvisPosX];
pelvisPosZ = [time,pelvisPosZ];
pelvisAngX = [time,pelvisAngX];
pelvisAngY = [time,pelvisAngY];
pelvisAngZ = [time,pelvisAngZ];
hipRAngX = [time,deg2rad(subj1.km.Rhip(:,1))];
hipRAngY = [time,deg2rad(subj1.km.Rhip(:,2))];
hipRAngZ = [time,-1*deg2rad(subj1.km.Rhip(:,3))]; 
hipLAngX = [time,deg2rad(subj1.km.Lhip(:,1))];
hipLAngY = [time,deg2rad(subj1.km.Lhip(:,2))];
hipLAngZ = [time,-1*deg2rad(subj1.km.Lhip(:,3))]; 

ankle_motionLY = [time,deg2rad(subj1.km.Lankle(:,2))];
ankle_motionLX = [time,deg2rad(subj1.km.Lankle(:,1))];
ankle_motionLZ = [time,-1*deg2rad(subj1.km.Lankle(:,3))];


ankle_motionR = [time,deg2rad(subj1.km.Rankle(:,2))];
knee_motionL = [time,deg2rad(subj1.km.Lknee(:,2))];
knee_motionR = [time,deg2rad(subj1.km.Rknee(:,2))];


% Static Motion Inputs
% % gaitPeriod = 2;
% gaitPeriod = 32.27;
% time = linspace(0,gaitPeriod,7)';
% time_mts = linspace(0,gaitPeriod,100);
% ankle_motion = deg2rad([0 0 0 0 0 0 0]');
% knee_motion = deg2rad([0, 0, 0, 0, 0, 0, 0]');
% hip_motion = deg2rad([0, 0, 0, 0, 0, 0, 0]');
% curveData = createSmoothTrajectory(ankle_motion,knee_motion,hip_motion,gaitPeriod);
% GRFv = linspace(0,700,7);

%% 


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

%% Open model
walkingRobot_LPVSFv6

%% Stiffness calculations
% order = 2;
% fc = 30; %hz cut off
% fs = length(tout)/.8;
% Wn = 2*fc/fs; %for forceplate data
% [b,a] = butter(order,Wn); % forceplate
% CF_filt = filtfilt(b,a,MTS_keel_contact.signals.values(:,7)); %symmetric lowpass filter of GRFS
% 
% 
% plot(MTS_position.signals.values,CF_filt);
% plot(MTS_position.time,MTS_position.signals.values);
% 
% plot(MTS_keel_contact.time,MTS_keel_contact.signals.values(:,7))
% plot(MTS_keel_contact.time(400:end),CF_filt(400:end));
% 
% keel_force = MTS_keel_contact.signals.values(:,7);
% mts_disp = ((MTS_position.signals.values - .042224)*1000);
% time = MTS_position.time;
