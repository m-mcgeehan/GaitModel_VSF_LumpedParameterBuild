%Startup script for lumped parameter VSF model

%% Clear the things
% clc
% clear
% close all

%% Add dependencies to the path

addpath(genpath('ModelingSimulation'), ...      % Modeling and simulation files
        genpath('Optimization'), ...            % Optimization files
        genpath('ControlDesign'), ...           % Control design files
        genpath('ReinforcementLearning'), ...   % Reinforcement learning files
        genpath('matlab scripts'), ...          %Directory for matlab scripts and functions
        genpath('Subject Bone Geometry'), ...   %Generic bone geometry files
        genpath('VSF CAD parts'),...            %VSF CAD files
        genpath('Experimental gait data'), ...  %Folder with experimental gait data
        genpath('SM_Contact_Forces_Lib'),...    %Contact forces library
        genpath('Libraries'));                  % Other dependencies

  
%% Model controls
GravityOn = 1; %1=on, 0=off...select zero for MTS simulation
fulcrum_position = 108; %mm
ModelControl = 1;%2 = v3d-defined kinematics, 1 = user-defined kinematics, 0 = static
frontalPlaneLock = 1;
coronalPlaneLock = 1;
LLside = 1; %1 = R, 0 = L

%% Load data
[fName,fPath] = uigetfile('.mat');
addpath(genpath(fPath));
load(fName);
name = extractBefore(fName,'.');
eval([('dynamic=') (convertStringsToChars(name))]); 
clear fName fPath name
%% Build skeletal model
ScaleBoneGeometryv3 %scale generic bone geometries to subject-specific anatomy based on marker data

%% Parameterize model
if isequal(LLside,1)
    GaitModelParams_VSF_R_2_24_20
else
    vsf_walker_paramsL
end

% edit vsf_walker_params_R_BodyParts3D; %run to edit model parameters
% edit vsf_walker_paramsL; %run to edit model parameters