%% Parameterize model
%Gravity
disp('Parameterizing model...');

if isequal(GravityOn,1) 
    gravity = -9.80665;
else
    gravity = 0;
end

%% Gait Model Geometry

%variable stiffness foot geometry
keel_length_total = 229; %mm
keel_thickness = 6.35; %mm
keel_width = 68.4; %mm
keel_length_rigid = 66; %mm
keel_length_tip = 78; %mm
keel_discretized_length = 10.5; % mm
keel_density = 6572.37; %kg/m^3

%General parameters
bone_color = [1 1 1];
density = 1000;
foot_density = 2000;
world_damping = 300; %0.25;
world_rot_damping = 1000; %0.25;
% if ~exist('actuatorType','var')
%     actuatorType = 1;
% % end

%ANATOMIC PARAMETERS

%Segment Length parameters
%Foot
foot_x = dynamic.static.segLength.footLX ;
foot_y = dynamic.static.segLength.footLY;
foot_z = dynamic.static.segLength.footLZ;

%Leg
leg_radius = 9.5;
shank_radius = 5.5;

lower_leg_length = dynamic.static.segLength.shankLengthL;
upper_leg_length = dynamic.static.segLength.legLengthL;
lower_leg_length_amputated = dynamic.static.segLength.shankLengthR;
res_shank_length = 100*lower_leg_length*.66; %residual leg length expressed as a % of intact shank length
pylon_radius = 1.5; %cm
socket2shank_offset = [-.5 -1 -1*(res_shank_length+5.5)]; %cm 2.5
pylon_length = 100*dynamic.static.segLength.shankLengthL + socket2shank_offset(3) + (dynamic.static.segLength.footLZ- 6.2); %cm
pylon_density = 5000; %kg/m^3

%Torso (cm)
torso_x = dynamic.static.segLength.pelvisX;
torso_y = dynamic.static.segLength.pelvisY;
torso_z = dynamic.static.segLength.pelvisZ;

%Mass parameters
SubjMass = dynamic.static.mass; %kg
mass_HAT = 0.628*SubjMass; % regression model limb masses from de Leva, 1996
mass_leg = .1416*SubjMass;
    volume_leg = ((pi/3)*(upper_leg_length))*((leg_radius/100)^2+(shank_radius/100)^2+(leg_radius/100*shank_radius/100));
    leg_density = mass_leg/volume_leg;
mass_shank = .0433*SubjMass;
mass_shank_amputated = mass_shank*((res_shank_length/lower_leg_length)/100);
    volume_shank = ((pi/3)*(lower_leg_length))*((shank_radius/100)^2+(shank_radius/100)^2+(shank_radius/100*shank_radius/100));
    shank_density = mass_shank/volume_shank;
mass_foot = .0137*SubjMass;
mass_total = mass_HAT+mass_leg*2+mass_shank+mass_shank_amputated+mass_foot+.8;%.8 = vsf mass
if (SubjMass-mass_total)>0
    mass_pylon = .5*(SubjMass-mass_total); %kg
    mass_socket = .5*(SubjMass-mass_total); %kg
else
    mass_pylon = .1;
    mass_socket = .1;
end

mass_total = mass_HAT+mass_leg*2+mass_socket+mass_pylon+mass_shank+mass_shank_amputated+mass_foot+.8;%.8 = vsf mass

%% Model assembly

% Joint offsets %socket to shank offset listed under segment length parameters
osf = 1; %Offset scale factor (used to adjust joint and segment offsets when scaling between subjects

world2pelvis = dynamic.static.offsets.World2Pelvis;
pelvis2HAT = dynamic.static.markers.ASISmidpt - dynamic.static.markers.pelvisOrigin;

%Right
% pelvis2HJCR = dynamic.static.offsets.Pelvis2HJCR;
legLength = 10*(norm(dynamic.static.markers.ASISmidpt - dynamic.static.markers.LMEP) + norm(dynamic.static.markers.LMEP - dynamic.static.markers.LHME));
pelvis2HJCR = (.1*[(11-(.063*legLength)); -1*(8+(.086*legLength)); (-9-(.078*legLength))]);
HJCR2KJCR = dynamic.static.offsets.FKneeRPos_LegR;
KJCR2shankR = dynamic.static.offsets.FKneeRPos_ShankR;

%Left
% pelvis2HJCL = dynamic.static.offsets.Pelvis2HJCL; 
pelvis2HJCL = (.1*[(11-(.063*legLength)); (8+(.086*legLength)); (-9-(.078*legLength))]);
HJCL2KJCL = dynamic.static.offsets.FKneeLPos_LegL;
KJCL2shankL = dynamic.static.offsets.FKneeLPos_ShankL;
shankL2AJCL = dynamic.static.offsets.FootLPos_ShankL;

foot_offset_vsf = [2.5 1.5 -3];
foot_offset = [-2.5 -1.5 4.5];
footL_anatomic_offset = [-1.8 0 1];

pelvisPos_init = [0 0 0];
ankleAngleR_init = [0 0 0];
ankleAngleL_init = [0 0 0];
hipAngleR_init = [0 0 0];
hipAngleL_init = [0 0 0];
kneeAngleR_init =[0 0 0];
kneeAngleL_init = [0 0 0];
pelvisAngle_init = [0 0 0];

pelvisOrigin2ASISR = dynamic.static.markers.RASI -dynamic.static.markers.pelvisOrigin;
pelvisOrigin2ASISL = dynamic.static.markers.LASI -dynamic.static.markers.pelvisOrigin;
pelvisOrigin2PSISR = dynamic.static.markers.RPSI -dynamic.static.markers.pelvisOrigin;
pelvisOrigin2PSISL = dynamic.static.markers.LPSI -dynamic.static.markers.pelvisOrigin;


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
motion_time_constant = .01; %(60/dynamic.framerate{1})/2; %0.025;

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
%     revolute_66_stiffness = k_inf; 
%     revolute_76_5_stiffness = k_inf;
%     revolute_87_stiffness = 6.75;
%     revolute_97_5_stiffness = k_inf;
%     revolute_108_stiffness = 6.75;
%     revolute_118_5_stiffness = k_inf;
%     revolute_129_stiffness = 6.75;
%     revolute_139_5_stiffness = k_inf;
%     revolute_150_stiffness = 6.75; 
%     revolute_160_5_stiffness = k_inf; 
%     revolute_171_stiffness = 6.75;
%     revolute_181_5_stiffness = k_inf;
%     revolute_192_stiffness = 6.75;
%     revolute_202_5_stiffness = k_inf;
%     revolute_213_stiffness = 6.75;
    revolute_66_stiffness = k_inf; 
    revolute_76_5_stiffness = krot;
    revolute_87_stiffness = k_inf;
    revolute_97_5_stiffness = krot;
    revolute_108_stiffness = k_inf;
    revolute_118_5_stiffness = krot;
    revolute_129_stiffness = k_inf;
    revolute_139_5_stiffness = krot;
    revolute_150_stiffness = k_inf; 
    revolute_160_5_stiffness = krot;
    revolute_171_stiffness = k_inf;
    revolute_181_5_stiffness = krot;
    revolute_192_stiffness = k_inf;
    revolute_202_5_stiffness =krot;
    revolute_213_stiffness = k_inf;
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
   revolute_damping = krot*.1;%5;
end

%% Motion Inputs

gaitPeriod = length(dynamic.km.Rankle(:,1))/dynamic.framerate{1};
time = linspace(0,gaitPeriod,length(dynamic.km.Rankle))';
if isequal(ModelControl,0)
        zeros = linspace(0,0,length(time))';
        
        pelvisPosX = [time,zeros];
        pelvisPosY = [time,zeros];  
        pelvisPosZ = [time,dynamic.km.PelvisCOMPos(:,3)];
        pelvisAngX = [time,zeros];
        pelvisAngY = [time,zeros];
        pelvisAngZ = [time,zeros];

        hipRAngX = [time,zeros];
        hipRAngY = [time,zeros];
        hipRAngZ = [time,zeros];
        hipLAngX = [time,zeros];
        hipLAngY = [time,zeros];
        hipLAngZ = [time,zeros];

        kneeRAngX = [time,zeros];
        kneeRAngY = [time,zeros];
        kneeRAngZ = [time,zeros];
        kneeLAngX = [time,zeros];
        kneeLAngY = [time,zeros];
        kneeLAngZ = [time,zeros];

        ankleLAngX = [time,zeros];
        ankleLAngY = [time,zeros];
        ankleLAngZ = [time,zeros];
        ankle_motionRX = -1*[time,deg2rad(dynamic.km.Rankle(:,1))];
            a = linspace(0,0,length(time));
        ankle_motionRY = [time,a'];  
        ankle_motionRZ = [time, -1*deg2rad(dynamic.km.Rankle(:,3))];
        timeDelay = 0;
        simTime = gaitPeriod;
    
    else
        if isequal(ModelControl,1)
        gaitPeriod = length(dynamic.km.Rankle(:,1))/dynamic.framerate{1};
        time = linspace(0,gaitPeriod,length(dynamic.km.Rankle))';

        pelvisPosX = [time,dynamic.km.PelvisCOMPos(:,1)];
        pelvisPosY = [time,dynamic.km.PelvisCOMPos(:,2)];  
        pelvisPosZ = [time,dynamic.km.PelvisCOMPos(:,3)];
        pelvisAngX = [dynamic.km.TS.pelvisAng_GCS.time, -1*deg2rad(dynamic.km.TS.pelvisAng_GCS.data(:,1))];
        pelvisAngY = [dynamic.km.TS.pelvisAng_GCS.time, -1*deg2rad(dynamic.km.TS.pelvisAng_GCS.data(:,2))];
        pelvisAngZ = [dynamic.km.TS.pelvisAng_GCS.time, deg2rad(dynamic.km.TS.pelvisAng_GCS.data(:,3))];

        hipRAngX = [dynamic.km.TS.RHipAngle.time, -1*deg2rad(dynamic.km.TS.RHipAngle.data(:,1))];
        hipRAngY = [dynamic.km.TS.RHipAngle.time, -1*deg2rad(dynamic.km.TS.RHipAngle.data(:,2))];
        hipRAngZ = [dynamic.km.TS.RHipAngle.time, deg2rad(dynamic.km.TS.RHipAngle.data(:,3))];
        hipLAngX = [dynamic.km.TS.LHipAngle.time, -1*deg2rad(dynamic.km.TS.LHipAngle.data(:,1))];
        hipLAngY = [dynamic.km.TS.LHipAngle.time, -1*deg2rad(dynamic.km.TS.LHipAngle.data(:,2))];
        hipLAngZ = [dynamic.km.TS.LHipAngle.time, deg2rad(dynamic.km.TS.LHipAngle.data(:,3))];

        kneeRAngX = [dynamic.km.TS.RKneeAngle.time, -1*deg2rad(dynamic.km.TS.RKneeAngle.data(:,1))];
        kneeRAngY = [dynamic.km.TS.RKneeAngle.time, -1*deg2rad(dynamic.km.TS.RKneeAngle.data(:,2))];
        kneeRAngZ = [dynamic.km.TS.RKneeAngle.time, deg2rad(dynamic.km.TS.RKneeAngle.data(:,3))];
        kneeLAngX = [dynamic.km.TS.LKneeAngle.time, -1*deg2rad(dynamic.km.TS.LKneeAngle.data(:,1))];
        kneeLAngY = [dynamic.km.TS.LKneeAngle.time, -1*deg2rad(dynamic.km.TS.LKneeAngle.data(:,2))];
        kneeLAngZ = [dynamic.km.TS.LKneeAngle.time, deg2rad(dynamic.km.TS.LKneeAngle.data(:,3))];

        ankleLAngX = [dynamic.km.TS.LAnkleAngle.time, -1*deg2rad(dynamic.km.TS.LAnkleAngle.data(:,1))];
        ankleLAngY = [dynamic.km.TS.LAnkleAngle.time, -1*deg2rad(dynamic.km.TS.LAnkleAngle.data(:,2))];
        ankleLAngZ = [dynamic.km.TS.LAnkleAngle.time, deg2rad(dynamic.km.TS.LAnkleAngle.data(:,3))];   
        ankle_motionRX = -1*[time,deg2rad(dynamic.km.Rankle(:,1))];
            a = linspace(0,0,length(time));
        ankle_motionRY = [time,a'];  
        ankle_motionRZ = [time, -1*deg2rad(dynamic.km.Rankle(:,3))];
        timeDelay = 0;
        simTime = gaitPeriod;

        else %V3d kinematics
            if isequal(ModelControl,2)
            gaitPeriod = length(dynamic.km.Rankle(:,1))/dynamic.framerate{1};
            time = linspace(0,gaitPeriod,length(dynamic.km.Rankle))';

            pelvisPosX = [time,dynamic.km.PelvisCOMPos(:,1)/100];
            pelvisPosY = [time,dynamic.km.PelvisCOMPos(:,2)/100];  
            pelvisPosZ = [time,dynamic.km.PelvisCOMPos(:,3)/100];
            pelvisAngX = [dynamic.km.v3d.TS.pelvisAng_GCS_v3d.time, deg2rad(dynamic.km.v3d.TS.pelvisAng_GCS_v3d.data(:,1))];
            pelvisAngY = [dynamic.km.v3d.TS.pelvisAng_GCS_v3d.time, deg2rad(dynamic.km.v3d.TS.pelvisAng_GCS_v3d.data(:,2))]; %-dynamic.static.km.v3d.pelvisAng_GCS_v3d(2)
            pelvisAngZ = [dynamic.km.v3d.TS.pelvisAng_GCS_v3d.time, deg2rad(dynamic.km.v3d.TS.pelvisAng_GCS_v3d.data(:,3))];

            hipRAngX = [dynamic.km.v3d.TS.RHipAngle.time, deg2rad(dynamic.km.v3d.TS.RHipAngle.data(:,1))];
            hipRAngY = [dynamic.km.v3d.TS.RHipAngle.time, deg2rad(dynamic.km.v3d.TS.RHipAngle.data(:,2))];
            hipRAngZ = [dynamic.km.v3d.TS.RHipAngle.time, deg2rad(dynamic.km.v3d.TS.RHipAngle.data(:,3))];
            hipLAngX = [dynamic.km.v3d.TS.LHipAngle.time, deg2rad(dynamic.km.v3d.TS.LHipAngle.data(:,1))];
            hipLAngY = [dynamic.km.v3d.TS.LHipAngle.time, deg2rad(dynamic.km.v3d.TS.LHipAngle.data(:,2))];
            hipLAngZ = [dynamic.km.v3d.TS.LHipAngle.time, deg2rad(dynamic.km.v3d.TS.LHipAngle.data(:,3))];

            kneeRAngX = [dynamic.km.v3d.TS.RKneeAngle.time, deg2rad(dynamic.km.v3d.TS.RKneeAngle.data(:,1))];
            kneeRAngY = [dynamic.km.v3d.TS.RKneeAngle.time, deg2rad(dynamic.km.v3d.TS.RKneeAngle.data(:,2))];
            kneeRAngZ = [dynamic.km.v3d.TS.RKneeAngle.time, deg2rad(dynamic.km.v3d.TS.RKneeAngle.data(:,3))];
            kneeLAngX = [dynamic.km.v3d.TS.LKneeAngle.time, deg2rad(dynamic.km.v3d.TS.LKneeAngle.data(:,1))];
            kneeLAngY = [dynamic.km.v3d.TS.LKneeAngle.time, deg2rad(dynamic.km.v3d.TS.LKneeAngle.data(:,2))];
            kneeLAngZ = [dynamic.km.v3d.TS.LKneeAngle.time, deg2rad(dynamic.km.v3d.TS.LKneeAngle.data(:,3))];

            ankleLAngX = [dynamic.km.v3d.TS.LAnkleAngle.time, deg2rad(dynamic.km.v3d.TS.LAnkleAngle.data(:,1))];
            ankleLAngY = [dynamic.km.v3d.TS.LAnkleAngle.time, deg2rad(dynamic.km.v3d.TS.LAnkleAngle.data(:,2)-dynamic.static.km.v3d.LAnkleAngle(2))];
            ankleLAngZ = [dynamic.km.v3d.TS.LAnkleAngle.time, deg2rad(dynamic.km.v3d.TS.LAnkleAngle.data(:,3))];   

            ankle_motionRX = -1*[time,deg2rad(dynamic.km.Rankle(:,1))];
                a = linspace(0,0,length(time));
            ankle_motionRY = [time,a'];  
            ankle_motionRZ = [time, -1*deg2rad(dynamic.km.Rankle(:,3))];

            timeDelay = 0;
            simTime = gaitPeriod; 
            else
                error('undefined model control');
            end
        end
end

if frontalPlaneLock == 1
%     pelvisAngX = [time,a'];
    hipRAngX = [time,a'];
    hipLAngX = [time,a'];
    kneeRAngX = [time,a'];
    kneeLAngX = [time,a'];
    ankleLAngX = [time,a'];
end
if coronalPlaneLock == 1
    pelvisPosY = [time,a'];  
%     pelvisAngZ = [time,a'];
    hipRAngZ = [time,a'];
    hipLAngZ = [time,a'];
    kneeRAngZ = [time,a'];
    kneeLAngZ = [time,a'];
    ankleLAngZ = [time,a'];
end
%% 

%% Contact model parameters

%% 
%Contact parameters
height_plane = 0.025;
plane_x = 25;
plane_y = 3;
mu_k = 0.3356;
mu_s = 0.1172;
mu_vth = 1.18;

%VSF

%VSF
contact_stiffness = 24120;
contact_damping = 15635;
contact_point_radius = .0011;
contact_color = [1 0 0];
contact_opc = 0;
heelCS_opc = 0.5;
ck_heel_VSF = 1034;
cd_heel_VSF = 3585.3;
xtrans = 12;
ztrans_heel_VSF = 2.99;
penDepth_heel_VSF = .0013;
penExp_heel_VSF = 7.66;

rad_heel_VSF = .0424;


%intact
contact_stiffness_intact_heel = 33145;
contact_damping_intact_heel = 49366;

contact_stiffness_intact_lateral = 5356.3;
contact_damping_intact_lateral = 4843.3;

contact_stiffness_intact_MTP = 3958.4;
contact_damping_intact_MTP = 74.36;

contact_stiffness_intact_forefoot = 13014;
contact_damping_intact_forefoot = 245.85;

contact_color_intact = [0 .5 .6];
contact_opc_intact = 0;
heel_intact_rad = .03;
lateral_intact_rad = .025;

pen_exp_heel_intact = 16.753;
pen_depth_heel_intact = 0.0203;
pen_exp_lateral_intact = 15.56;
pen_depth_lateral_intact = 0.0059;
% init_height = foot_z + dynamic.static.segLength.legLengthL + dynamic.static.segLength.shankLengthL + ...
%     dynamic.static.segLength.pelvisZ + height_plane/2;
% init_height = foot_z + dynamic.static.segLength.legLengthL + dynamic.static.segLength.shankLengthL + ...
%     dynamic.static.segLength.pelvisZ + height_plane/2;
init_height = dynamic.km.PelvisCOMPos(1,3)*100;
init_pos = [dynamic.km.PelvisCOMPos(1,1)*100 dynamic.km.PelvisCOMPos(1,2)*100 dynamic.km.PelvisCOMPos(1,3)*100];

world_model_damping = 1;
world_model_stiffness = 1;

bone_color = [1 1 1];

limb_geom_opc = 0.2;

marker_color = [0.0 0.4 1.0];
marker_radius = 0.5; %cm

%% 
pyrD = 1272;
pyrK = 2192.6;
socketD_rot = 4980;
socketD_trans = 9936.5;
socketK_rot = 9945;
socketK_trans = 19915;

VSF_heel_prisK = 32546;
VSF_heel_prisD = 32245;
heel_prisK = 11092;
heel_prisD = 18743;

%%
vsf_fc = 45;
intact_fc = 45;
%% Set up model geometry file paths
boneGeomPath = fileparts(which('pelvis_scaled.stl'));
fpPelvis = fullfile(boneGeomPath,'pelvis_scaled.stl');
fpHAT = fullfile(boneGeomPath,'HAT_scaled');
fpFemurR = fullfile(boneGeomPath,'FemurR_scaled');
fpFemurL = fullfile(boneGeomPath,'FemurL_scaled');
fpShankR = fullfile(boneGeomPath,'ShankR_cut_scaled.stl');
fpShankL = fullfile(boneGeomPath,'ShankL_scaled.stl');
fpFootL_hindfoot = fullfile(boneGeomPath,'FootL_hindfoot_scaled.stl');
fpFootL_forefoot = fullfile(boneGeomPath,'FootL_forefoot_scaled.stl');

%% Additional stuffs *Clean this up later
% contact_d_heel_mid = 1500;
% contact_k_heel_mid = 2000;
% contact_k_heel_ant = 2000;
% contact_d_heel_ant = 1000;
% contact_k_heel_post = 5000;
% contact_d_heel_post = 4500;
% contact_k_heel_pc = 10000;
% contact_d_heel_pc = 8000;



% from LHC optimzation
% contact_k_heel_post = 520;
% contact_d_heel_post = 1748;
% contact_k_heel_pc = 448;
% contact_d_heel_pc = 720;
% 
% contact_r_heel = .001;

% heel_rad = 0.03;
% heel_pen = .01;
% pen_exp = 1;

% From intact optimization
% contact_damping_intact_MTP = 118.7;
% contact_damping_intact_forefoot = 160.9;
% contact_damping_intact_heel = 2910.9;
% contact_damping_intact_lateral = 4036.5;

% contact_stiffness_intact_MTP = 118.7;
% contact_stiffness_intact_MTP = 10118.7;

% contact_stiffness_intact_forefoot = 160.9;
% contact_stiffness_intact_forefoot = 25600.9;
% 
% contact_stiffness_intact_heel = 2910.9;
% contact_stiffness_intact_lateral = 5874.4;
% 
% pen_depth_heel_intact = .00751;
% pen_depth_lateral_intact = .00574;
% pen_exp_heel_intact = 1.4297;
% pen_exp_lateral_intact = 1.3741;
% heel_intact_rad = .03;
% lateral_intact_rad = .03;
% 
% For HCM2 heel contact model
% ck_heel_VSF = 50000;
% cd_heel_VSF = 500;

vars_GaitModel3_11_20

Torque = Simulink.Variant('actuationMode == 1');
Motion = Simulink.Variant('actuationMode == 0');

disp('Parameterizing model...Done');

%% Open model
disp('Loading GaitModelVSF_4_1_20');

GaitModelVSF_4_1_20;

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
