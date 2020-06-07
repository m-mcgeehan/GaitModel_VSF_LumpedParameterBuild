function penalty = simGaitModel(params,mdlName,scaleFactor,period,~)
% Cost function for robot walking optimization
% Copyright 2017 The MathWorks, Inc.

    % Load parameters into function workspace
    vsf_walker_params_skeletonR
    % Apply variable scaling
    params = scaleFactor*params;
    
    % Extract simulation inputs from parameters
    N = numel(params)/3;
    
%     contact_stiffness = 4000;
%     contact_damping =.6*contact_stiffness;


    % Simulate the model
    simout = sim(mdlName,'SrcWorkspace','current','FastRestart','on');          

    % Unpack logged data
    
    
    wAvg = mean(simout.yout{1}.Values.Data);
    xEnd = simout.yout{2}.Values.Data(end);
    tEnd = simout.tout(end);

    % Calculate penalty from logged data
    %   Distance traveled without falling 
    %   (ending simulation early) increases reward
    positiveReward = sign(xEnd)*xEnd^2 * tEnd;
    %   Angular velocity and trajectory aggressiveness 
    %   (number of times the derivative flips signs) decreases reward
    %   NOTE: Set lower limits to prevent divisions by zero
    
    
order = 2;
fs = length(measR.normal_force.time)/gaitPeriod; 
fc = 20;%50 hz cut off
Wn_fp = 2*fc/fs; 
[bfp,afp] = butter(order,Wn_fp); 
simForceR = filtfilt(bfp,afp,measR.normal_force.data); %symmetric lowpass filter of normal force
timeSim = linspace(0,gaitPeriod,length(measR.normal_force.time));
time1 = linspace(0,gaitPeriod,length(subj1.forces.FP1(:,3)));

fSim = simForceR(5726:9083);
fExp = interp1(time1,subj1.forces.FP2(:,3),timeSim)';
fExp = fExp(5726:9083);


deltaSignal = abs(fSim - fExp);
percentageDifference = deltaSignal ./ fExp; % Percent by element.
negativeReward = mean(percentageDifference); % Average percentage over all elements.
    
   
    penalty = -positiveReward/negativeReward;        
    
end

