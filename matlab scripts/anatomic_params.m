contact_stiffness = 100000;
density = 1000;
HAT_density = 8000;

contact_damping = 5000;
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

