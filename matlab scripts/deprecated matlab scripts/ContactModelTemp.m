contact_stiffness = 4000;
contact_damping =.6*contact_stiffness;

contact_stiffness_heel = 100;
contact_damping_heel =.99*contact_stiffness;

contact_stiffness_intact = 3000;
contact_damping_intact =.6*contact_stiffness;


mu_k = 0.9;
mu_s = 0.99;
mu_vth = 2;
height_plane = 0.025;
plane_x = 25;
plane_y = 3;
contact_point_radius = .01;
contact_color = [0 0 0];
contact_opc = 0;


contact_color_intact = [1 0 0];
contact_opc_intact = 0;

init_height = foot_z + lower_leg_length + upper_leg_length + ...
              torso_z/2 + torso_offset_z + height_plane/2 - 5.4;


world_model_damping = 1;
world_model_stiffness = 1;

bone_color = [1 1 1];