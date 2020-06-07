function costFcn = forceCostFunction(measR, subj1)

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
costFcn = mean(percentageDifference); % Average percentage over all elements.


% plot(fSim); hold on
% plot(fExp);

