% Main script for optimizing the gait of a walking robot model
% Copyright 2017 The MathWorks, Inc.

%% Set initial parameters

% Choose model name
mdlName = 'walkingRobot_LPVSF_R_skeleton_v2'; % Main model
% mdlName = 'RobotWalkerModified'; % CAD imported model

% Flags to speed up simulation
accelFlag = true;
parallelFlag = false;

% Joint actuator type for optimization
% 1 = motion | 2 = torque | 3 = motor
actuatorType = 1;

% To reduce the search space, scale the angle waypoints and solve the
% optimization algorithm with integer parameters
scalingFactor = 1;        % Scaling from degrees to integer parameters

% Uncomment to create initial conditions from scratch
numPoints = length(measR.normal_force.data); % Number of GRF data points 
gaitPeriod = simTime;           % Period of walking gait [s]
p0 = zeros(1,numPoints*3);  % Create zero motion initial conditions

% Uncomment to load initial gait from previous optimization results
% startFile = 'optimizedData_17Jul17_1021';
% [p0,gaitPeriod,numPoints] = createInitialConditions(startFile,scalingFactor);

%% Set optimization options
opts = optimoptions('ga');
opts.Display = 'iter';
opts.MaxGenerations = 10;
opts.PopulationSize = 10;
opts.InitialPopulationMatrix = repmat(p0,[5 1]); % Add copies of initial gait
opts.PlotFcn = @gaplotbestf; % Add progress plot of fitness function
opts.UseParallel = parallelFlag;

%% Set bounds and constraints
% Upper and lower angle bounds
upperBnd = [10000*ones(1,numPoints), ... % stiffness limits
            9999*ones(1,numPoints)] ... % damping limits
            /scalingFactor;
lowerBnd = [1000*ones(1,numPoints), ... % stiffness limits
            100*ones(1,numPoints)] ... % damping limits
            /scalingFactor;

%% Run commands to set up parallel/accelerated simulation
doSpeedupTasks;

%% Run optimization
costFcn = @(p)simulateGaitModel(p,mdlName,scalingFactor,gaitPeriod,actuatorType);
disp(['Running optimization. Population: ' num2str(opts.PopulationSize) ...
      ', Max Generations: ' num2str(opts.MaxGenerations)])
[pFinal,reward] = ga(costFcn,numPoints*3,[],[],[],[], ... 
                     lowerBnd,upperBnd,[],1:numPoints*3,opts);
disp(['Final reward function value: ' num2str(-reward)])

%% Save results to MAT-file
% Convert from optimization integer search space to trajectories in radians
pScaled = scalingFactor*pFinal;
time = linspace(0,gaitPeriod,numPoints+1)';
% ankle_motion = deg2rad([pScaled(1:numPoints) pScaled(1)]');
% knee_motion = deg2rad([pScaled(numPoints+1:2*numPoints) pScaled(numPoints+1)]');
% hip_motion = deg2rad([pScaled(2*numPoints+1:3*numPoints) pScaled(2*numPoints+1)]');
% curveData = createSmoothTrajectory(ankle_motion,knee_motion,hip_motion,gaitPeriod);
% outFileName = ['optimizedData_' datestr(now,'ddmmmyy_HHMM')];
% save(outFileName,'curveData','reward','gaitPeriod','time', ... 
%                  'hip_motion','knee_motion','ankle_motion');
% plotSmoothTrajectory;

%% Cleanup
bdclose(mdlName);
if parallelFlag
   delete(gcp('nocreate')); 
end