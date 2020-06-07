%% Calculate CoP from force magnitude and sphere coordinates

%% Process individual contact sphere normal forces
%% open and process experimental CoP data
%Align coordinate systems (simscape wrt v3d: x = y, y = -x, z = z)
x = [0 1 0;...
    -1 0 0;...
    0 0 1];

COP_VSF = COP_VSF{1}';
    indices = find(abs(COP_VSF)>1);
    COP_VSF(indices) = 0;
    
for i = 1:length(COP_VSF(1,:))  %Rotate COP dat and convert to cm
    COP_VSF(:,i) = (x*(COP_VSF(:,i))*100);
end

% a = max(COP_VSF(3,:));
COP_VSFap = (COP_VSF(3,:))'; clear a



%%

fs = 200; %FS of downsampled forces = 200Hz
buffer = fs*.25; %number of samples to pad on either side of stance onset/offset
time = linspace(0,(length(dynamic.forces.FP2)/fs),length(dynamic.forces.FP2));
phaseShift = .00; %phase shift (s) from analog low-pass filter at 40 Hz in simulation

idx1 = find(dynamic.forces.FP2(:,3) > 1,1,'first');
idx2 = find(dynamic.forces.FP2(:,3)> 1,1,'last');

% VSFGRFR_exp = timeseries((sqrt(dynamic.forces.FP2(:,3).^2 + dynamic.forces.FP2(:,2).^2 + ...
%     dynamic.forces.FP2(:,1).^2)),time+phaseShift);

VSFGRFR_exp = timeseries((sqrt(dynamic.forces.FP2(idx1-buffer:idx2+buffer,3).^2 + dynamic.forces.FP2(idx1-buffer:idx2+buffer,2).^2 + ...
    dynamic.forces.FP2(idx1-buffer:idx2+buffer,1).^2)),(time(idx1-buffer:idx2+buffer))+phaseShift);
VSFCOP_exp = timeseries(abs(COP_VSFap(idx1-buffer:idx2+buffer))*10,(time(idx1-buffer:idx2+buffer))+phaseShift);


VSFFn_interp = timeseries(interp1(SDOSimTest_Log.tout, SDOSimTest_Log.VSF_Fn_unfiltered.data,VSFGRFR_exp.time, 'nearest'),VSFGRFR_exp.time); %Interpolate sim force
VSFFf_interp = timeseries(interp1(SDOSimTest_Log.tout, SDOSimTest_Log.VSF_Ff_unfiltered.data,VSFGRFR_exp.time, 'nearest'),VSFGRFR_exp.time); %Interpolate sim force
VSFGRFR_interp = timeseries((VSFFn_interp.data+VSFFf_interp.data),VSFGRFR_exp.time);

% plot(VSFGRFR_exp); hold on; plot(VSFGRFR_interp);

idx1 = find(VSFGRFR_exp.data > 50,1,'first');
idx2 = find(VSFGRFR_exp.data > 50,1,'last');
% plot(VSFGRFR_exp.data(idx1-1:idx2+1));
VSFCOP_exp = timeseries(VSFCOP_exp.data(idx1-1:idx2+1),VSFCOP_exp.time(idx1-1:idx2+1));


idx1 = find(VSFGRFR_interp.data > 50,1,'first');
idx2 = find(VSFGRFR_interp.data > 50,1,'last');
if idx1<2
    idx1 = 2;
else
    
end

% plot(VSFGRFR_interp.data(idx1-1:idx2+1));
VSFGRFR_interp = timeseries(VSFGRFR_interp.data(idx1-1:idx2+1),VSFGRFR_interp.time(idx1-1:idx2+1));


forceSignals_Fn = {'heelCenter_Fn' 'heelPostM_Fn' 'heelPostL_Fn' 'x66M_Fn' 'x66L_Fn' 'x76_Fn' 'x129_Fn' 'x139M_Fn' 'x139L_Fn' 'x150M_Fn' 'x150L_Fn'...
    'x160M_Fn' 'x160L_Fn' 'x171M_Fn' 'x171L_Fn' 'x181M_Fn' 'x181L_Fn' 'x192M_Fn' 'x192L_Fn' 'x202M_Fn' 'x202L_Fn' 'x213M_Fn' 'x213L_Fn' 'x229_Fn'};

forceSignals_Ff = {'heelCenter_Ff' 'heelPostM_Ff' 'heelPostL_Ff' 'x66M_Ff' 'x66L_Ff' 'x76_Ff' 'x129_Ff' 'x139M_Ff' 'x139L_Ff' 'x150M_Ff' 'x150L_Ff'...
    'x160M_Ff' 'x160L_Ff' 'x171M_Ff' 'x171L_Ff' 'x181M_Ff' 'x181L_Ff' 'x192M_Ff' 'x192L_Ff' 'x202M_Ff' 'x202L_Ff' 'x213M_Ff' 'x213L_Ff' 'x229_Ff'};

heelCenter_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_24_.data;
heelPostM_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_19_.data;
heelPostL_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_20_.data;
x66M_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_21_.data;
x66L_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_22_.data;
x76_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_23_.data;
x129_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_10_.data;
x139M_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_13_.data;
x139L_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_11_.data;
x150M_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_12_.data;
x150L_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_14_.data;
x160M_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_16_.data;
x160L_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_15_.data;
x171M_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_18_.data;
x171L_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_17_.data;
x181M_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_1_.data;
x181L_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_5_.data;
x192M_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_2_.data;
x192L_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_3_.data;
x202M_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_4_.data;
x202L_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_6_.data;
x213M_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_8_.data;
x213L_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_7_.data;
x229_Fn = SDOSimTest_Log.VSF_normForces.Normal_Force__signal_9_.data;

heelCenter_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_24_.data;
heelPostM_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_19_.data;
heelPostL_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_20_.data;
x66M_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_21_.data;
x66L_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_22_.data;
x76_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_23_.data;
x129_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_10_.data;
x139M_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_13_.data;
x139L_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_11_.data;
x150M_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_12_.data;
x150L_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_14_.data;
x160M_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_16_.data;
x160L_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_15_.data;
x171M_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_18_.data;
x171L_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_17_.data;
x181M_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_1_.data;
x181L_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_5_.data;
x192M_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_2_.data;
x192L_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_3_.data;
x202M_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_4_.data;
x202L_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_6_.data;
x213M_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_8_.data;
x213L_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_7_.data;
x229_Ff = SDOSimTest_Log.VSF_frictForces.Friction_Force__signal_9_.data;

order = 2;
fs = 200;
fc = 20;%40 hz low pass cut off
Wn_fp = 2*fc/fs; 
[b,a] = butter(order,Wn_fp); 

% force_total = filtfilt(b,a,VSFGRFR_interp.data);
% plot(force_total);

for i = 1:length(forceSignals_Fn)
    signal_Fn = forceSignals_Fn{i};
    eval(['signalData_Fn =',signal_Fn,';']);
    signalInterp_Fn = interp1(SDOSimTest_Log.tout, signalData_Fn,VSFGRFR_interp.time, 'nearest');
    signalFilt_Fn = filtfilt(b,a,signalInterp_Fn);
    normForce.(signal_Fn) = signalFilt_Fn;
    
    signal_Ff = forceSignals_Ff{i};
    eval(['signalData_Ff =',signal_Ff,';']);
    signalInterp_Ff = interp1(SDOSimTest_Log.tout, signalData_Ff,VSFGRFR_interp.time, 'nearest');
    signalFilt_Ff = filtfilt(b,a,signalInterp_Ff);
    frictForce.(signal_Ff) = signalFilt_Ff;
    
    extractAfter(string(FILE_NAME{1}),"VSF\"); 
    name = strcat((extractBefore(signal_Ff,'_')),'_total');
    totalForce.(name) = signalFilt_Ff + signalFilt_Fn;
end


%% Define x coordinates for each sphere
heelCenterX = 3;
heelPostX = 21;
x66X = 71.2;
x76X = 82.2;
x129X = 130;
x139X = 144.7;
x150X = 155.7;
x160X = 166.2;
x171X = 176.7;
x181X = 187.2;
x192X = 197.7;
x202X = 208.2;
x213X = 218.7;
x229X = 229;

%% Calculate CoP as weighted sum of product of Fn and x coordinate of each sphere

order = 2;
fs = 200;
fc = 10;%40 hz low pass cut off
Wn_fp = 2*fc/fs; 
[b,a] = butter(order,Wn_fp);


% test = ((heelCenterX.*normForce.heelCenter)+(heelPostX.*normForce.heelPostL))./force_total;

% copX = filtfilt(b,a,(heelCenterX*normForce.heelCenter)+(heelPostX*normForce.heelPostM)+(heelPostX*normForce.heelPostL)+(x66X*normForce.x66L)+(x66X*normForce.x66M)...
%     +(x76X*normForce.x76)+(x129X*normForce.x129)+(x139X*normForce.x139L)+(x139X*normForce.x139M)...
%     +(x150X*normForce.x150L)+(x150X*normForce.x150M)+(x160X*normForce.x160L)+(x160X*normForce.x160M)+(x171X*normForce.x171L)+(x171X*normForce.x171M)...
%     +(x181X*normForce.x181L)+(x181X*normForce.x181M)+(x192X*normForce.x192L)+(x192X*normForce.x192M)+...
%     (x202X*normForce.x202L)+(x202X*normForce.x202M)+(x213X*normForce.x213L)+(x213X*normForce.x213M)+(x229X*normForce.x229)/force_total);


% copX = filtfilt(b,a,((heelCenterX.*normForce.heelCenter)+(heelPostX.*normForce.heelPostM)+(heelPostX.*normForce.heelPostL)+(x66X.*normForce.x66L)+(x66X.*normForce.x66M)...
%     +(x76X.*normForce.x76)+(x129X.*normForce.x129)+(x139X.*normForce.x139L)+(x139X.*normForce.x139M)...
%     +(x150X.*normForce.x150L)+(x150X.*normForce.x150M)+(x160X.*normForce.x160L)+(x160X.*normForce.x160M)+(x171X.*normForce.x171L)+(x171X.*normForce.x171M)...
%     +(x181X.*normForce.x181L)+(x181X.*normForce.x181M)+(x192X.*normForce.x192L)+(x192X.*normForce.x192M)+...
%     (x202X.*normForce.x202L)+(x202X.*normForce.x202M)+(x213X.*normForce.x213L)+(x213X.*normForce.x213M)+(x229X.*normForce.x229))./force_total);

FnSum = (normForce.heelCenter_Fn+normForce.heelPostM_Fn...
    +normForce.heelPostL_Fn+normForce.x66L_Fn+normForce.x66M_Fn+normForce.x76_Fn+normForce.x129_Fn+normForce.x139L_Fn+normForce.x139M_Fn...
    +normForce.x150L_Fn+normForce.x150M_Fn+normForce.x160L_Fn+normForce.x160M_Fn+normForce.x171L_Fn+normForce.x171M_Fn...
    +normForce.x181L_Fn+normForce.x181M_Fn+normForce.x192L_Fn+normForce.x192M_Fn+...
    normForce.x202L_Fn+normForce.x202M_Fn+normForce.x213L_Fn+normForce.x213M_Fn+normForce.x229_Fn);

FfSum = (frictForce.heelCenter_Ff+frictForce.heelPostM_Ff...
    +frictForce.heelPostL_Ff+frictForce.x66L_Ff+frictForce.x66M_Ff+frictForce.x76_Ff+frictForce.x129_Ff+frictForce.x139L_Ff+frictForce.x139M_Ff...
    +frictForce.x150L_Ff+frictForce.x150M_Ff+frictForce.x160L_Ff+frictForce.x160M_Ff+frictForce.x171L_Ff+frictForce.x171M_Ff...
    +frictForce.x181L_Ff+frictForce.x181M_Ff+frictForce.x192L_Ff+frictForce.x192M_Ff+...
    frictForce.x202L_Ff+frictForce.x202M_Ff+frictForce.x213L_Ff+frictForce.x213M_Ff+frictForce.x229_Ff);

FSum = FnSum + FfSum;

% copX = timeseries(filtfilt(b,a,((heelCenterX.*normForce.heelCenter)+(heelPostX.*normForce.heelPostM)+(heelPostX.*normForce.heelPostL)+(x66X.*normForce.x66L)+(x66X.*normForce.x66M)...
%     +(x76X.*normForce.x76)+(x129X.*normForce.x129)+(x139X.*normForce.x139L)+(x139X.*normForce.x139M)...
%     +(x150X.*normForce.x150L)+(x150X.*normForce.x150M)+(x160X.*normForce.x160L)+(x160X.*normForce.x160M)+(x171X.*normForce.x171L)+(x171X.*normForce.x171M)...
%     +(x181X.*normForce.x181L)+(x181X.*normForce.x181M)+(x192X.*normForce.x192L)+(x192X.*normForce.x192M)+...
%     (x202X.*normForce.x202L)+(x202X.*normForce.x202M)+(x213X.*normForce.x213L)+(x213X.*normForce.x213M)+(x229X.*normForce.x229))./FnSum),(time(idx1-buffer:idx2+buffer)));

copX = timeseries(filtfilt(b,a,((heelCenterX.*totalForce.heelCenter_total)+(heelPostX.*totalForce.heelPostM_total)+(heelPostX.*totalForce.heelPostL_total)+(x66X.*totalForce.x66L_total)+(x66X.*totalForce.x66M_total)...
    +(x76X.*totalForce.x76_total)+(x129X.*totalForce.x129_total)+(x139X.*totalForce.x139L_total)+(x139X.*totalForce.x139M_total)...
    +(x150X.*totalForce.x150L_total)+(x150X.*totalForce.x150M_total)+(x160X.*totalForce.x160L_total)+(x160X.*totalForce.x160M_total)+(x171X.*totalForce.x171L_total)+(x171X.*totalForce.x171M_total)...
    +(x181X.*totalForce.x181L_total)+(x181X.*totalForce.x181M_total)+(x192X.*totalForce.x192L_total)+(x192X.*totalForce.x192M_total)+...
    (x202X.*totalForce.x202L_total)+(x202X.*totalForce.x202M_total)+(x213X.*totalForce.x213L_total)+(x213X.*totalForce.x213M_total)+(x229X.*totalForce.x229_total))./FSum),VSFGRFR_interp.time);

copX101 = timeseries(imresize(copX.data,[101 1]),imresize(copX.time,[101 1]));
VSFCOP_exp101 = timeseries(imresize(VSFCOP_exp.data,[101 1]),imresize(VSFCOP_exp.time,[101 1]));

copX_RMSE = sqrt(immse(copX101.data,VSFCOP_exp101.Data));
copX_r2 = rsquared(copX101.data,VSFCOP_exp101.Data);


%%

plot(copX101.data); hold on;
plot(VSFCOP_exp101.data);


% plot(copX(58:178)/10); hold on;
% plot(abs(COP_VSFap(311:431)));
% a = find(VSFCOP_exp.data > 1,1,'last');
% plot(copX.time(1:a),copX.data(1:a)); hold on;
% plot(VSFCOP_exp.time(1:a),abs(VSFCOP_exp.data(1:a)));

% 
% plot(copX(1:a)/10); hold on;
% plot(VSFCOP_exp.time,abs(VSFCOP_exp.data));
% 
% 
% COP_RMSE = sqrt(immse((copX(58:178)/10),2),abs(COP_VSFap(311:431)));
% COP_r1 = rsquared(copX(58:178)/10),abs(COP_VSFap(311:431));


