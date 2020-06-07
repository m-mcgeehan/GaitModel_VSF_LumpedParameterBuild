% Walking Robot Parameters
% Copyright 2017 The MathWorks, Inc.

%% General parameters
density = 1000;
foot_density = 2000;
world_damping = 0.25;
world_rot_damping = 0.25;
if ~exist('actuatorType','var')
    actuatorType = 1;
end

%% Inputs
gaitPeriod = 0.8;
time = linspace(0,gaitPeriod,7)';
ankle_motion = deg2rad([-7.5 10 10 5 0 -10 -7.5]');
knee_motion = deg2rad([10, -5, 2.5, -10, -10, 15, 10]');
hip_motion = deg2rad([-10, -7.5, -15, 10, 15, 10, -10]');
curveData = createSmoothTrajectory(ankle_motion,knee_motion,hip_motion,gaitPeriod);

%% Contact/friction parameters
contact_stiffness = 2500;
contact_damping = 100;
mu_k = 0.6;
mu_s = 0.8;
mu_vth = 0.1;
height_plane = 0.025;
plane_x = 25;
plane_y = 3;
contact_point_radius = 1e-4;

%% Foot parameters
foot_x = 8;
foot_y = 6;
foot_z = 1;
foot_offset = [-1 0 0];

%% Leg parameters
leg_radius = 0.75;
lower_leg_length = 10;
upper_leg_length = 10;

%% Torso parameters
torso_y = 10;
torso_x = 5;
torso_z = 8;
torso_offset_z = -2;
torso_offset_x = -0.5;
init_height = foot_z + lower_leg_length + upper_leg_length + ...
              torso_z/2 + torso_offset_z + height_plane/2;

%% Joint parameters
joint_damping = 1;
joint_stiffness = 1;
motion_time_constant = 0.01; %0.025;

%% Joint controller parameters
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

%% Electric motor parameters
motor_resistance = 1;
motor_constant = 0.02;
motor_inertia = 0;
motor_damping = 0;
motor_inductance = 1.2e-6;
gear_ratio = 50;

%% added from LPVSF model

gravity_on = 1; %1=on, 0=off...select zero for MTS simulation
fulcrum_position = 66; %mm

%Gravity
if isequal(gravity_on,1) 
    gravity = -9.80665;
else
    gravity = 0;
end

%Pylon geometry
pylon_radius = 15; %mm
pylon_length = 200; %mm
pylon_density = 50000; %kg/m^3

%variable stiffness foot geometry
keel_length_total = 229; %mm
keel_thickness = 6.35; %mm
keel_width = 68.4; %mm
keel_length_rigid = 66; %mm
keel_length_tip = 78; %mm
keel_discretized_length = 10.5; %(keel_length_total-(keel_length_tip+keel_length_rigid))/keel_descretization; %mm
keel_density = 6572.37; %kg/m^3

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

foot_offset_vsf = [-2.7 31 -3.3];

lower_leg_length_amputated = 10;

mts_OuterRadius = .01; %m
mts_length = .068; %m
mts_density = .001; %kg/m^3

%MTS kinematics
%for motion control (comment out for force control):
max_disp = 26.894; %mm
mts_displacement = linspace(42.224,42.224+max_disp,100); %[5;20;35;50;65;80;95]; %mm
mts_contact_damping = 100; %N/m/s
mts_contact_stiffness = 1000000; %N/m
%for force control:
mts_force = linspace(0,800,100);

%% remove this
contact_stiffness = 100000;
density = 1000;
contact_damping = 5000
mu_k = 0.9;
mu_s = 0.9;
mu_vth = 0.01;
height_plane = 0.025;
plane_x = 25;
plane_y = 3;
contact_point_radius = 1e-4;


foot_x = 22;
foot_y = 7;
foot_z = 7;
foot_offset_vsf = [-2.7 31 -3.3];
% foot_offset = [-1 0 0];
foot_offset = [-7 0 0];

%Pylon geometry
pylon_radius = 15; %mm
pylon_length = 200; %mm
pylon_density = 50000; %kg/m^3

leg_radius = 0.75;
lower_leg_length = 10;
upper_leg_length = 10;



%Leg parameters (cm)

leg_radius = 6;
lower_leg_length = 40;
upper_leg_length = 40;
lower_leg_length_amputated = 10;

%Torso parameters (cm)
torso_y = 27;
torso_x = 16;
torso_z = 16;
torso_offset_z = -2;
torso_offset_x = -0.5;
init_height = foot_z + lower_leg_length + upper_leg_length + ...
              torso_z/2 + torso_offset_z + height_plane/2;
init_height = init_height - 3;


