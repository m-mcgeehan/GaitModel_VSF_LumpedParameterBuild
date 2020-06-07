VSFGRFR_sim = SDOSimTest_Log.VSF_CF_resultant;
intactGRFR_sim = SDOSimTest_Log.intact_CF;

VSFGRFR_interp = interp1(VSFGRFR_sim.time, VSFGRFR_sim.data,VSFGRFR_expShift.time, 'spline', 'extrap'); 
VSFGRFR_interp = VSFGRFR_interp(1:90);
VSFGRFR_expShift = timeseries(VSFGRFR_expShift.data(1:90),VSFGRFR_expShift.time(1:90));

VSF_RMSE = sqrt(immse(VSFGRFR_expShift.data,VSFGRFR_interp));
VSFGRFR_r2 = rsquared(VSFGRFR_interp,VSFGRFR_expShift.data,2);

intactGRFR_interp = interp1(intactGRFR_sim.time, intactGRFR_sim.data, intactGRFR_expShift.time,'spline','extrap');
intact_RMSE = sqrt(immse(intactGRFR_expShift.data,intactGRFR_interp));
intactGRFR_r2 = rsquared(intactGRFR_interp,intactGRFR_expShift.data,2);

%Impulse
VSFGRFR_exp_impulse = cumtrapz(VSFGRFR_expShift.time,VSFGRFR_expShift.data);
VSFGRFR_sim_impulse = cumtrapz(VSFGRFR_expShift.time,VSFGRFR_interp);
VSF_impulse_RMSE = sqrt(immse(VSFGRFR_exp_impulse,VSFGRFR_sim_impulse));
VSFGRFR_impulse_r2 = rsquared(VSFGRFR_exp_impulse,VSFGRFR_sim_impulse,2);

intact_exp_impulse = cumtrapz(intactGRFR_expShift.time,intactGRFR_expShift.data);
intact_sim_impulse = cumtrapz(intactGRFR_expShift.time,intactGRFR_interp);
intact_impulse_RMSE = sqrt(immse(intact_exp_impulse,intact_sim_impulse));
intact_impulse_r2 = rsquared(intact_exp_impulse,intact_sim_impulse,2);

figure;
subplot(2,2,1); plot(VSFGRFR_expShift.time, VSFGRFR_interp,'k'); hold on; plot(VSFGRFR_expShift,'r');...
    legend('simulation','experimental'); title(['VSF GRFR R^2 = ', num2str(VSFGRFR_r2), ' RMSE =', num2str(VSF_RMSE), 'N']);
subplot(2,2,2); plot(intactGRFR_expShift.time, intactGRFR_interp,'k'); hold on; plot(intactGRFR_expShift,'r'); legend('simulation','experimental');
    legend('simulation','experimental'); title(['Intact GRFR R^2 = ', num2str(intactGRFR_r2), ' RMSE =', num2str(intact_RMSE), 'N']);
subplot(2,2,3); plot(VSFGRFR_expShift.time, VSFGRFR_sim_impulse,'k'); hold on; plot(VSFGRFR_expShift.time,VSFGRFR_exp_impulse,'r');...
    legend('simulation','experimental'); title(['VSF Impulse R^2 = ', num2str(VSFGRFR_impulse_r2), ' RMSE =', num2str(VSF_impulse_RMSE), 'N/m/s']);
subplot(2,2,4); plot(intactGRFR_expShift.time, intact_sim_impulse,'k'); hold on; plot(intactGRFR_expShift.time,intact_exp_impulse,'r');...
    legend('simulation','experimental'); title(['Intact Impulse R^2 = ', num2str(intact_impulse_r2), ' RMSE =', num2str(intact_impulse_RMSE), 'N/m/s']);
sgtitle(strcat('Resultant GRF for: ', dynamic.trial));